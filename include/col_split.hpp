#include "common/common.hpp"
#include <cassert>
#include <common.hpp>
#include <FL_table.hpp>
#include <cstddef>
#include <fstream>
#include <r_index.hpp>

#include <iostream>
#include <sdsl/io.hpp>
#include <sdsl/int_vector.hpp>
#include <malloc_count.h>

#ifndef _COL_SPLIT_TABLE_HH
#define _COL_SPLIT_TABLE_HH

namespace Options {
    /* How to split bwt runs based on cols */
    enum class Mode {
        Default, // Split all
        Tunneled, // Split only until FL_mapping diverges
        RunAligned // Split only if cols are run aligned
        // TunneledAll // Split only until FL_mapping diverges, but keep keep stepping in case it converges
    };

    /* How to handle overlapping runs */
    enum class Overlap {
        Append,  // Keep existing col_run, split to include non-overlapping parts of col run to be added
        Split,  // Split on overlaps (keep all non-overlapping parts even if it segments into multiple runs)
        Remove, // Remove the entirety of any col run with overlaps
        Ignore // Remove run which overlaps with existing col run (FCFS)
        // Largest -> maximize coverage
    };
}

template <class FL_t = FL_table>
class col_split {
public:
    using Mode = Options::Mode;
    using Overlap = Options::Overlap;

    col_split() {}

    col_split(FL_t& t) : tbl(t), n(tbl.size()), r(tbl.runs()) {}

    col_split(FL_t& t, Mode m, Overlap o = Overlap::Append) : tbl(t), n(tbl.size()), r(tbl.runs()), mode(m), overlap(o) {}

    col_split(FL_t& t, Overlap o, Mode m = Mode::Default) : tbl(t), n(tbl.size()), r(tbl.runs()), mode(m), overlap(o) {}

    void set_mode(Mode m) {
        mode = m;
    }

    void set_overlap(Overlap o) {
        overlap = o;
    }

    // Split bases on Co-linear positions and lengths; assuming they have uniform height N
    // Split rate: split at this rate while forward stepping to carve out runs (all run aligned cols are split anyway)
    void split(vector<ulint> &col_len, vector<ulint> &col_pos, ulint N, int split_rate = 1) {
        assert(col_len.size() == col_pos.size());

        status("Initializing vectors");
        bitvec col_run_bv(n);
        // TODO use O(r') data structure for col_ids (linked list)
        vector<ulint> col_ids(n, 0);

        // Mark all ids which would have been split if split rate was 1
        // used to idenitfy run-aligned cols
        phantom_ids = phantom_ids_t(n, split_rate);
        // TODO destroy phantom_ids on exit
        // TODO find run aligned during computation to avoid O(run_len) scans
        // TODO for memory efficiency, use phantom_id bools and flip (keep col_ids true, use phantom_ids for false 0s, etc)
        status();

        status("FL Stepping");
        size_t run_start = 0;
        size_t c_id = 2; // 0 denotes to id, 1 denotes overlap
        vector<ulint>::iterator c_pos = col_pos.begin();
        vector<ulint>::iterator c_len = col_len.begin();

        // first pass finds col run boundaries
        for (size_t i = 0; i < tbl.runs(); ++i) {
            ulint run_len = tbl.get_length(i);

            // TODO Fix - N on mumemto and col_bwt side
            while (c_pos != col_pos.end() && *c_pos - N < run_start + run_len) {
                range rg = {i, *c_pos - N - run_start, N};
                vector<range> FL_ranges = FL_range(rg);

                bool skip_non_tunnel = (mode == Mode::Tunneled) && FL_ranges.size() > 1;
                for (size_t j = 0; j < *c_len && !skip_non_tunnel; ++j) {
                    vector<range> next_ranges;
                    for (size_t k = 0; k < FL_ranges.size(); ++k) {
                        rg = FL_ranges[k];
                        ulint idx = tbl.get_idx(rg.interval);
                        ulint col_span_start = idx + rg.offset;
                        ulint col_span_end = idx + rg.offset + rg.height;

                        // Has boundary overlap, so skip computation for ignore/remove modes
                        bool skip_overlap = (overlap == Overlap::Remove || overlap == Overlap::Ignore) && 
                            (col_ids[col_span_start] > 1 || col_ids[col_span_end - 1] > 1);
                        if (skip_overlap)
                        {
                            overlap_skip(col_ids, col_run_bv, col_span_start, col_span_end);
                        }
                        else {
                            ulint last_id = (col_span_start > 0) ? col_ids[col_span_start - 1] : 0;
                            for (size_t l = col_span_start; l < col_span_end && !skip_overlap; ++l) {
                                // force split if overlapping
                                bool force_split = false;
                                
                                if (col_ids[l] >= 1) {
                                    switch (overlap) {
                                        case Overlap::Split:
                                            col_ids[l] = 1; // mark overlaps
                                            phantom_ids.set(l, 1);
                                            col_run_bv.unset(l); // Erase any col runs here; they have overlaps

                                            /*
                                            * (col_ids[l] == 1) --> overlapping, need to split regardless of split rate
                                            * (last_id == 1) --> last run was overlapping, need to split regardless of split rate
                                            */
                                            force_split = col_ids[l] == 1 || last_id == 1;
                                            break;
                                        case Overlap::Remove: 
                                        case Overlap::Ignore:
                                            // TODO Can check the boundaries first to save time
                                            overlap_skip(col_ids, col_run_bv, col_span_start, col_span_end, l);
                                            skip_overlap = true;
                                            break;
                                        default: // Overlap::Append does nothing, leave unchanged
                                            break;
                                    }
                                }
                                else {
                                    if (j % split_rate == 0) {
                                        col_ids[l] = c_id;
                                    }
                                    else {
                                        col_ids[l] = 0;
                                    }

                                    phantom_ids.set(l, c_id);
                                }

                                // Skip the rest of the loop if this col run is no longer being considered
                                if (!skip_overlap) {
                                    bool split_rate_met = j % split_rate == 0 && (col_ids[l] == c_id || last_id == c_id);
                                    if (last_id != col_ids[l] && (split_rate_met || force_split)) {
                                        col_run_bv.set(l);
                                    }
                                    else if (last_id == col_ids[l]) {
                                        col_run_bv.unset(l);
                                    }

                                    last_id = col_ids[l];
                                }
                            }
                            
                            // Close run if needed
                            if (!skip_overlap && col_span_end < n) {
                                last_id = col_ids[col_span_end - 1];
                                bool id_mismatch = (last_id != col_ids[col_span_end]);
                                bool split_rate_met = (j % split_rate == 0) && (last_id == c_id);
                                bool need_run_end = false;
                                if (overlap == Overlap::Split) {
                                    // End the current run if:
                                    // i) we are at the end of a run of overlapping characters and the next is not overlapping
                                    // ii) we are at the end of a run of non-overlapping characters and the next id does not match
                                    bool force_split = (last_id > 1 || col_ids[col_span_end] > 0);
                                    need_run_end = id_mismatch && (force_split || split_rate_met);
                                }
                                else {
                                    need_run_end = id_mismatch && split_rate_met;
                                }

                                if (need_run_end) {
                                    col_run_bv.set(col_span_end);
                                }
                                // If we can continue the run, do so
                                else if (last_id == col_ids[col_span_end]) {
                                    col_run_bv.unset(col_span_end);
                                }
                            }
                        }

                        vector<range> curr_ranges = FL_range(rg);
                        next_ranges.insert(next_ranges.end(), curr_ranges.begin(), curr_ranges.end());
                    }
                    FL_ranges = next_ranges;
                    skip_non_tunnel = (mode == Mode::Tunneled) && FL_ranges.size() > 1;
                }
                ++c_pos;
                ++c_len;
                ++c_id;
            }
            run_start += run_len;
        }
        status();

        // second pass resolves col run boundaries with bwt run boundaries
        find_col_runs(col_ids, col_run_bv, split_rate);
    }

    size_t save_full(string filename) {
        size_t written_bytes = 0;

        written_bytes += serialize_col_runs(col_runs, filename);

        std::ofstream out_ids(filename + ".col_ids");
        for (size_t i = 0; i < col_run_ids.size(); ++i) {
            out_ids.write(reinterpret_cast<const char*>(&col_run_ids[i]), RW_BYTES);
            written_bytes += RW_BYTES;
        }
        return written_bytes;
    }

    /* Space saving option to mod the values before saving */
    size_t save(string filename, int id_bits=ID_BITS) {
        assert (id_bits <= sizeof(ulint) * 8);
        ulint id_max = 1ULL << id_bits;
        ulint id_bytes = (id_bits + 7) / 8; // bytes needed to store id

        size_t written_bytes = 0;

        written_bytes += serialize_col_runs(col_runs, filename);

        std::ofstream out_ids(filename + ".col_ids");
        for (size_t i = 0; i < col_run_ids.size(); ++i) {
            size_t id = col_run_ids[i] % id_max;
            out_ids.write(reinterpret_cast<const char*>(&id), id_bytes);
            written_bytes += id_bytes;
        }
        return written_bytes;
    }

private:
    typedef struct {
        ulint interval;
        ulint offset;
        ulint height;
    } range;

    struct bitvec {
        bit_vector bv;
        size_t set_count;

        bitvec() : bv(), set_count(0) {}

        bitvec(size_t size) : bv(size), set_count(0) {}

        bool operator[](size_t idx) const { return bv[idx]; }
        bit_vector get_bv() const { return bv; }

        size_t size() const { return bv.size(); }
        size_t set_bits() const { return set_count; }
        size_t unset_bits() const { return size() - set_bits(); }

        void set(size_t idx) {
            set_count += !bv[idx];
            bv[idx] = 1;
        }

        void unset(size_t idx) {
            set_count -= bv[idx];
            bv[idx] = 0;
        }
    };

    struct phantom_ids_t{
        vector<ulint> p_ids;
        bool split_rate_used = false;

        phantom_ids_t() : p_ids(), split_rate_used(false) {}

        phantom_ids_t(size_t n, int split_rate) : split_rate_used(split_rate > 1) {
            if (split_rate_used) {
                p_ids = vector<ulint>(n, 0);
            }
        }

        vector<ulint>& get_vec() {
            if (!split_rate_used) {
                error("Attempting to access p_ids when split_rate_used is false");
            }
            return p_ids;
        }

        void set(size_t idx, ulint id) {
            if (split_rate_used) {
                p_ids[idx] = id;
            }
        }
    };

    FL_t& tbl;
    ulint n = 0;
    ulint r = 0;
    Mode mode = Mode::Default;
    Overlap overlap = Overlap::Append;

    sd_vector<> col_runs;
    vector<ulint> col_run_ids;
    phantom_ids_t phantom_ids;

    vector<range> FL_range(range rg) {
        vector<range> FL_dest_ranges;

        // FL steps may diverge to different positions even within the same run
        while (rg.height > 0) {
            const auto& [FL_interval, FL_offset] = tbl.FL(rg.interval, rg.offset);

            if (rg.offset + rg.height > tbl.get_length(rg.interval))  {
                ulint covered = tbl.get_length(rg.interval) - rg.offset;
                FL_dest_ranges.push_back({FL_interval, FL_offset, covered});
                rg.height -= covered;
                rg.offset = 0;
            }
            else {
                FL_dest_ranges.push_back({FL_interval, FL_offset, rg.height});
                rg.height = 0;
            }
            // Move to the next run; the FL mapping is broken into multiple runs
            ++rg.interval;
        }
        return FL_dest_ranges;
    }

    // Overload to skip the current position (meaning still at col_span_start)
    void overlap_skip(vector<ulint> &col_ids, bitvec &col_run_bv, size_t col_span_start, size_t col_span_end) {
        overlap_skip(col_ids, col_run_bv, col_span_start, col_span_end, col_span_start);
    }

    void overlap_skip(vector<ulint> &col_ids, bitvec &col_run_bv, size_t col_span_start, size_t col_span_end, size_t col_span_curr) {
        switch (overlap) {
            case Overlap::Remove:
                overlap_remove(col_ids, col_run_bv, col_span_start, col_span_end, col_span_curr);
                break;
            case Overlap::Ignore:
                overlap_ignore(col_ids, col_run_bv, col_span_start, col_span_end, col_span_curr);
                break;
            default:
                break;
        }
    }

    void overlap_remove(vector<ulint> &col_ids, bitvec &col_run_bv, size_t col_span_start, size_t col_span_end, size_t col_span_curr) {
        // Remove overlapping run moving upwards (only occurs if col_span_start is new overlap)
        if (col_ids[col_span_start] > 1) {
            size_t curr_pos = col_span_start;
            size_t colliding_id = col_ids[curr_pos];
            while (curr_pos > 0 && col_ids[curr_pos - 1] == colliding_id) {
                col_ids[curr_pos - 1] = 1;
                phantom_ids.set(curr_pos - 1, 1);
                col_run_bv.unset(curr_pos - 1);
                --curr_pos;
            }
            if (curr_pos > 0 && col_ids[curr_pos - 1] > 1) {
                if (col_ids[curr_pos - 1] > 1) {
                    col_run_bv.set(curr_pos);
                }
                else {
                    col_run_bv.unset(curr_pos);
                }
            }
        }

        // Remove overlapping run moving downwards (only occurs if col_span_end - 1 is new overlap)      
        if (col_ids[col_span_end - 1] > 1) {
            size_t curr_pos = col_span_end;
            size_t colliding_id = col_ids[curr_pos];
            while (curr_pos < n && col_ids[curr_pos] == colliding_id) {
                col_ids[curr_pos] = 1;
                phantom_ids.set(curr_pos, 1);
                col_run_bv.unset(curr_pos);
                ++curr_pos;
            }
            if (curr_pos < n && col_ids[curr_pos + 1] <= 1) {
                col_run_bv.unset(curr_pos);
            }
        }

        // Mark current run as overlapping
        for (size_t curr_pos = col_span_start; curr_pos < col_span_end; ++curr_pos) {
            col_ids[curr_pos] = 1;
            phantom_ids.set(curr_pos, 1);
            col_run_bv.unset(curr_pos);
        }
        if (col_span_start > 0 && col_ids[col_span_start - 1] > 1) {
            col_run_bv.set(col_span_start);
        }
    }

    void overlap_ignore(vector<ulint> &col_ids, bitvec &col_run_bv, size_t col_span_start, size_t col_span_end, size_t col_span_curr) {
        // Ignore current col run by undoing work so far
        for (size_t curr_pos = col_span_start; curr_pos < col_span_curr; ++curr_pos) {
            col_ids[curr_pos] = 0;
            phantom_ids.set(curr_pos, 0);
            col_run_bv.unset(curr_pos);
        }
        if (col_span_start > 0 && col_ids[col_span_start - 1] > 1) {
            col_run_bv.set(col_span_start);
        }
    }

    // Finds col runs and computes col_run_ids and col_runs
    void find_col_runs(vector<ulint> &col_ids, bitvec &col_run_bv, int split_rate) {
        #ifdef PRINT_STATS
        size_t col_chars = 0;
        size_t num_col_bits = 0;
        size_t last_idx = 0;
        #endif

        sd_select col_select;
        // Don't need bitvector, only aligns to existing runs
        // TODO can borrow existing sparse bv from FL
        if (mode == Mode::RunAligned) {
            col_run_bv = bitvec(n);
        }
        // with these settings, we can compress the bitvector only once; print status now
        else if (split_rate == 1) {
            status("Creating runs bitvector");
            col_runs = sd_vector(col_run_bv.get_bv());
            col_select = sd_select(&col_runs);
            status();
        }
        else {
            col_runs = sd_vector(col_run_bv.get_bv());
            col_select = sd_select(&col_runs);
        }

        status("Creating IDs vector");
        col_run_ids = vector<ulint>();
        size_t set_bits = col_run_bv.set_bits();

        // To keep track of L_runs
        const vector<ulint> &run_ids = (split_rate > 1) ? phantom_ids.get_vec() : col_ids;
        size_t curr_L_run = 0;
        size_t curr_L_len = tbl.get_L_length(curr_L_run);
        size_t curr_L_pos = 0;
        ulint last_id = 0;

        // TODO refactor to not add runs if not needed
        auto process_bwt_runs = [&](size_t idx) {
            while (curr_L_run < tbl.runs() && curr_L_pos < idx) {
                if (mode == Mode::RunAligned || !col_run_bv[curr_L_pos]) {
                    bool run_aligned = false;
                    bool intersects_col_split = (idx < n) && (idx < curr_L_pos + curr_L_len);
                    // only check if run aligned mode or no id to assign, and ensure we do not intersect a colored split
                    bool check_run_aligned = mode == Mode::RunAligned || (last_id == 0 && !intersects_col_split);
                    if (check_run_aligned && run_ids[curr_L_pos] > 1 && (run_ids[curr_L_pos] == run_ids[curr_L_pos + curr_L_len - 1])) {
                        ulint curr_id = run_ids[curr_L_pos];
                        run_aligned = true;
                        for (size_t j = curr_L_pos + 1; j < curr_L_pos + curr_L_len - 1 && run_aligned; ++j) {
                            run_aligned = (run_ids[j] == curr_id);
                        }
                    }

                    if (run_aligned) {
                        col_run_bv.set(curr_L_pos);
                        add_col_run_id(run_ids[curr_L_pos]);
                        
                        #ifdef PRINT_STATS
                        col_chars += curr_L_len;
                        ++num_col_bits;
                        #endif
                    }
                    else {
                        ulint curr_id = (mode == Mode::RunAligned) ? 0 : last_id;
                        col_run_bv.set(curr_L_pos);
                        add_col_run_id(curr_id);
                    }
                }

                ++curr_L_run;
                if (curr_L_run < tbl.runs()) {
                    curr_L_pos = tbl.get_L_pos(curr_L_run);
                    curr_L_len = tbl.get_L_length(curr_L_run);
                }
            }
        };

        if (mode != Mode::RunAligned) {
            for (size_t i = 1; i <= set_bits; ++i) {
                size_t idx = col_select(i);

                // if subsampling, look for run aligned runs to add (does not add runs compared BWT runs)
                if (split_rate > 1) {
                    process_bwt_runs(idx);
                }

                // Add other non-run aligned
                add_col_run_id(col_ids[idx]);
                last_id = col_ids[idx];

                #ifdef PRINT_STATS
                ++num_col_bits;
                if (i > 1 && col_ids[last_idx] > 1) {
                    col_chars += idx - last_idx;
                }
                last_idx = idx;
                #endif
            }
            #ifdef PRINT_STATS
            if (col_ids[last_idx] > 1) {
                col_chars += n - last_idx;
                ++num_col_bits;
            }
            #endif
        }
        if (split_rate > 1 || mode == Mode::RunAligned) {
            process_bwt_runs(n);
        }
        status();

        // Need to remake runs bv if changes were made
        if (split_rate > 1 || mode == Options::Mode::RunAligned) {
            status("Creating runs bitvector");
            col_runs = sd_vector(col_run_bv.get_bv());
            status();
        }

        #ifdef PRINT_STATS
        cout << "Col ids: " << num_col_bits << std::endl;
        cout << "Col chars: " << col_chars << std::endl;
        #endif
    }

    void add_col_run_id(ulint id) {
        col_run_ids.push_back((id > 1) ? id - 1 : 0);
    }

    size_t serialize_col_runs(sd_vector<> &sdv, std::string filename) {
        size_t written_bytes = 0;

        std::ofstream out_runs(filename + ".col_runs");
        written_bytes += col_runs.serialize(out_runs);
        out_runs.close();

        return written_bytes;
    }
};

#endif /* end of include guard: _COL_SPLIT_HH */