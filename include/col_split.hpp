#include "common/common.hpp"
#include <cassert>
#include <common.hpp>
#include <FL_table.hpp>
#include <cstddef>
#include <fstream>
#include <r_index.hpp>

#include <sdsl/io.hpp>
#include <sdsl/int_vector.hpp>
#include <malloc_count.h>

#ifndef _COL_SPLIT_TABLE_HH
#define _COL_SPLIT_TABLE_HH

namespace Options {
    /* How to split bwt runs based on cols */
    enum class Mode {
        All, // Split all
        Tunneled, // Split only until FL_mapping diverges
    };

    /* How to handle overlapping runs */
    enum class Overlap {
        Append,  // Keep existing col_run, split to include non-overlapping parts of col run to be added
        Truncate  // Stop at overlap
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

    col_split(FL_t& t, Overlap o, Mode m = Mode::All) : tbl(t), n(tbl.size()), r(tbl.runs()), mode(m), overlap(o) {}

    void set_mode(Mode m) {
        mode = m;
    }

    void set_overlap(Overlap o) {
        overlap = o;
    }

    // Split bases on Co-linear positions and lengths; assuming they have uniform height N
    // Split rate: split at this rate while forward stepping to carve out runs (all run aligned cols are split anyway)
    void split(vector<ulint> &col_len, vector<ulint> &col_pos, ulint N, int split_rate = 1, int gap_rate = 20) {
        assert(col_len.size() == col_pos.size());

        vector<ulint> ids_per_start;
        if (gap_rate > 0) {
            ids_per_start = find_mum_blocks(col_len, col_pos, N, gap_rate);
        }

        status("Initializing n sized bitvector");
        bitvec mark_start_bv(n);
        status();

        auto collect_boundaries = [&] (size_t col_span_start, size_t _c_id, size_t _height) {
            mark_start_bv.set(col_span_start);
        };

        auto FL_loop = [&] (auto func) {
            size_t run_start = 0;
            size_t curr_col = 0;
            size_t c_id = 2; // 0 denotes to id, 1 denotes overlap
            if (gap_rate > 0) {
                c_id = ids_per_start[curr_col];
            }
            vector<ulint>::iterator c_pos = col_pos.begin();
            vector<ulint>::iterator c_len = col_len.begin();
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

                            if (j % split_rate == 0) 
                            {
                                ulint idx = tbl.get_idx(rg.interval);
                                ulint col_span_start = idx + rg.offset;
                                // Do some work here depending on the pass
                                func(col_span_start, c_id, rg.height);
                            }
                            
                            vector<range> curr_ranges = FL_range(rg);
                            next_ranges.insert(next_ranges.end(), curr_ranges.begin(), curr_ranges.end());
                        }
                        FL_ranges = next_ranges;
                        skip_non_tunnel = (mode == Mode::Tunneled) && FL_ranges.size() > 1;
                    }
                    ++c_pos;
                    ++c_len;
                    ++curr_col;
                    if (gap_rate > 0) {
                        c_id = ids_per_start[curr_col];
                    }
                    else {
                        ++c_id;
                    }
                }
                run_start += run_len;
            }
        };

        status("FL Stepping to find interval boundaries");
        FL_loop(collect_boundaries);
        status();

        status("Creating new vectors for Col-IDs");
        vector<std::pair<ulint,ulint>> marked_ids = vector<std::pair<ulint,ulint>>(mark_start_bv.set_bits(), {0,0});
        bit_rank rank_start(&mark_start_bv.get_bv());
        status();

        auto collect_ids = [&] (size_t col_span_start, size_t c_id, size_t height) {
            // For default mode, keep the largest height if overlap
            if (mode == Mode::All && mark_start_bv[col_span_start]) {
                size_t curr_rank = rank_start(col_span_start);
                auto [existing_id, existing_height] = marked_ids[curr_rank];
                size_t max_height = std::max(existing_height, height);
                size_t max_id = (existing_height >= height) ? existing_id : c_id;
                marked_ids[rank_start(col_span_start)] = {max_id, max_height};
            }
            // Tunnel case is easy, just mark the col run
            else {
                marked_ids[rank_start(col_span_start)] = {c_id, height};
            }
        };

        // Consider collecting the col starting positions to speed up second pass
        status("FL Stepping to fill Col-ID Values");
        FL_loop(collect_ids);
        status();

        // final pass resolves col run boundaries for overlaps and bwt run boundaries
        find_col_runs(mark_start_bv, marked_ids);
    }

    // size_t save(string filename) {
    //     size_t written_bytes = 0;

    //     written_bytes += serialize_col_runs(col_runs, filename);

    //     std::ofstream out_ids(filename + ".col_ids");
    //     written_bytes += write_vec(col_run_ids, out_ids);
    //     return written_bytes;
    // }

    /* Space saving option to mod the values before saving */
    size_t save(string filename, int id_bits = ID_BITS, bool sparse = true) {
        assert (id_bits <= sizeof(ulint) * 8);
        ulint id_max = bit_max(id_bits);
        ulint id_bytes = bits_to_bytes(id_bits); // bytes needed to store id

        size_t written_bytes = 0;

        written_bytes += serialize_col_runs(col_runs, filename, sparse);

        std::ofstream out_ids(filename + ".col_ids");
        for (size_t i = 0; i < col_run_ids.size(); ++i) {
            size_t id = col_run_ids[i];
            if (id >= id_max) {
                id = (id % (id_max - 1)) + 1;
            }
            out_ids.write(reinterpret_cast<const char*>(&id), id_bytes);
            written_bytes += id_bytes;
        }
        return written_bytes;
    }

    // size_t save_ids(string filename, int id_bits) {
    //     assert (id_bits <= sizeof(ulint) * 8);
    //     ulint id_max = bit_max(id_bits);
    //     ulint id_bytes = bits_to_bytes(id_bits); // bytes needed to store id

    //     size_t written_bytes = 0;

    //     std::ofstream out_ids(filename + ".col_ids");
    //     for (size_t i = 0; i < col_run_ids.size(); ++i) {
    //         size_t id = col_run_ids[i];
    //         if (id >= id_max) {
    //             id = (id % (id_max - 1)) + 1;
    //         }
    //         out_ids.write(reinterpret_cast<const char*>(&id), id_bytes);
    //         written_bytes += id_bytes;
    //     }
    //     return written_bytes;
    // }

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
        bit_vector& get_bv() { return bv; }

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

    FL_t& tbl;
    ulint n = 0;
    ulint r = 0;
    Mode mode = Mode::All;
    Overlap overlap = Overlap::Append;

    bit_vector col_runs;
    vector<ulint> col_run_ids;

    vector<ulint> find_mum_blocks(const vector<size_t> &col_len, const vector<size_t> &col_pos, ulint N, ulint gap_rate) {
        status("Assigning Col-IDs by finding co-linear MUM blocks");
        bitvec sa_starts_bv(n);
        vector<ulint> ids_per_start = vector<ulint>(col_pos.size());

        vector<ulint>::const_iterator c_pos = col_pos.begin();
        while (c_pos != col_pos.end()) {
            sa_starts_bv.set(*c_pos - N);
            ++c_pos;
        }
    
        bit_rank first_rank = bit_rank(&sa_starts_bv.get_bv());
        bit_select first_select = bit_select(&sa_starts_bv.get_bv());

        std::pair<ulint, ulint> curr_pos = {0, 0};
        size_t c_id = 2;
        ids_per_start[0] = c_id++; // fix weird error
        size_t visited = 0;
        size_t curr_mum_len = 0; 
        size_t curr_gap = 0;
        bool in_mum = 0;
        #ifdef PRINT_STATS
        size_t steps = 0;
        #endif
        while (visited <= col_pos.size()) {
            ulint idx = tbl.get_idx(curr_pos.first) + curr_pos.second;
            ulint last_rank = first_rank(idx+1);
            ulint last_idx = (last_rank > 0) ? first_select(last_rank) : 0;
            bool see_mum = (last_rank > 0 && idx - last_idx < N);

            if (see_mum) {
                if (!in_mum && curr_gap > gap_rate) {
                    ++c_id;
                }

                ids_per_start[last_rank] = c_id;
                curr_mum_len = col_len[last_rank];
                curr_gap = 0;
                ++visited;
                in_mum = true;
            }
            else if (in_mum) {
                if (curr_mum_len > 1) {
                    --curr_mum_len;
                } else {
                    in_mum = false;
                    curr_mum_len = 0;
                    curr_gap = 0;
                }
            }
            else {
                ++curr_gap;
            }

            curr_pos = tbl.FL(curr_pos.first, curr_pos.second);
            #ifdef PRINT_STATS
            ++steps;
            #endif
        }
        status();

        #ifdef PRINT_STATS
        cout << "IDs used: " << c_id - 1 << endl;
        cout << "  # Cols: " << col_pos.size() << endl;
        cout << " Visited: " << visited - 1 << endl;
        cout << "   Steps: " << steps << endl;
        cout << "     n/d: " << n/N << endl;
        #endif

        return ids_per_start;
    }

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

    void add_col_run_id(ulint id) {
        col_run_ids.push_back((id > 1) ? id - 1 : 0);
    }

    // Finds col runs and computes col_run_ids and col_runs
    void find_col_runs(bitvec &idx_start_bv, vector<std::pair<ulint,ulint>> &idx_col_ids) {
        if (idx_start_bv.set_bits() == 0) {
            return;
        }

        status("Loading BWT run support");
        col_runs = bit_vector(n, 0);
        bit_select start_select(&idx_start_bv.get_bv());
        bit_rank start_rank(&idx_start_bv.get_bv());

        sdsl::sd_vector<> bwt_runs = tbl.get_L_heads();
        sd_select bwt_run_select(&bwt_runs);
        status();

        struct interval {
            ulint start;
            ulint end;
            ulint id;

            bool operator<(const interval &other) const {
                return (end < other.end) || (end == other.end && start < other.start);
            }

            bool operator>(const interval &other) const {
                return (end > other.end) || (end == other.end && start > other.start);
            }
        };

        std::priority_queue<interval, vector<interval>, greater<interval>> open_intervals;
        size_t run_cursor = 1;
        ulint curr_bwt_pos = bwt_run_select(run_cursor);
        ulint last_id = 0;

        auto update_bwt_pos = [&](size_t idx, ulint id = 0) {
            while (run_cursor <= tbl.runs() && curr_bwt_pos < idx) {
                col_runs[curr_bwt_pos] = 1;
                add_col_run_id(last_id);
                ++run_cursor;
                curr_bwt_pos = bwt_run_select(run_cursor);
            }
            if (curr_bwt_pos == idx) {
                ++run_cursor;
                curr_bwt_pos = bwt_run_select(run_cursor);
            }
            last_id = id;
        };

        auto update_col_ranges = [&](size_t idx) {
            while (!open_intervals.empty() && open_intervals.top().end <= idx) {
                interval e = open_intervals.top();
                open_intervals.pop();
                if (open_intervals.size() == 1 && open_intervals.top().end > e.end) {
                    update_bwt_pos(e.end, open_intervals.top().id);
                    col_runs[e.end] = 1;
                    add_col_run_id(open_intervals.top().id);
                }
                else if (open_intervals.empty() && e.end < idx) {
                    update_bwt_pos(e.end, 0);
                    col_runs[e.end] = 1;
                    add_col_run_id(0);
                }
            }
        };

        status("Finding intervals of col sub-runs");
        for (size_t i = 1; i <= idx_start_bv.set_bits(); ++i) {
            size_t curr_start = start_select(i);
            auto [curr_col_id, curr_col_height] = idx_col_ids[start_rank(curr_start)];

            update_col_ranges(curr_start);

            open_intervals.push({curr_start, curr_start + curr_col_height, curr_col_id});
            if (open_intervals.size() == 1 && curr_col_id > 1) {
                update_bwt_pos(curr_start, curr_col_id);
                col_runs[curr_start] = 1;
                add_col_run_id(curr_col_id);
            }
        }
        update_col_ranges(n);
        update_bwt_pos(n, 0);
        status();

        #ifdef PRINT_STATS
        status("Finding stats on procedure");
        size_t col_chars = 0;
        size_t col_id_runs = 0;
        
        bit_rank col_rank(&col_runs);
        bit_select col_select(&col_runs);

        size_t set_bits = col_rank(n);
    
        size_t last_idx = col_select(1);
        last_id = col_run_ids[0];
        for (size_t i = 1; i <= set_bits; ++i) {
            size_t curr_idx = col_select(i);
            size_t curr_id = col_run_ids[i - 1];
            if (last_id >= 1) {
                ++col_id_runs;
                col_chars += curr_idx - last_idx;
            }
            last_idx = curr_idx;
            last_id = curr_id;
        }
        if (last_id >= 1) {
            col_chars += n - last_idx;
            ++col_id_runs;
        }
        status();

        cout << "Col runs: " << col_id_runs << std::endl;
        cout << "Total runs: " << set_bits << std::endl;
        cout << "Col chars: " << col_chars << std::endl;
        #endif
    }

    size_t serialize_col_runs(bit_vector &bv, std::string filename, bool sparse = true) {
        size_t written_bytes = 0;

        if (sparse) {
            sdsl::sd_vector<> sd_bv(bv);

            std::ofstream out_runs(filename + ".col_runs");
            written_bytes += sd_bv.serialize(out_runs);
            out_runs.close();
        } else {
            std::ofstream out_runs(filename + ".col_runs.bv");
            written_bytes += bv.serialize(out_runs);
            out_runs.close();
        }

        return written_bytes;
    }
};

#endif /* end of include guard: _COL_SPLIT_HH */