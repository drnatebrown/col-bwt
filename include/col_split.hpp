#include "common/common.hpp"
#include <cassert>
#include <common.hpp>
#include <FL_table.hpp>
#include <cstddef>
#include <fstream>
#include <optional>
#include <r_index.hpp>
#include <unordered_map>

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
    };

    /* How to handle overlapping runs */
    enum class Overlap {
        Append,  // Keep existing col_run, split to include non-overlapping parts of col run to be added
        Split  // Split on overlaps (keep all non-overlapping parts even if it segments into multiple runs)
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
        bitvec col_run_bv(n, N);
        std::unordered_map<ulint, ulint> col_ids;
        status();

        auto set_id = [&col_run_bv, &col_ids](ulint idx, ulint id) {
            col_run_bv.set(idx);
            col_ids[idx] = id;
        };

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

                        if (j % split_rate == 0) 
                        {
                            ulint idx = tbl.get_idx(rg.interval);
                            ulint col_span_start = idx + rg.offset;
                            ulint col_span_end = idx + rg.offset + rg.height;

                            ulint mark_id = c_id;
                            bool mark = true;

                            auto pred = col_run_bv.predecessor(col_span_start + 1);
                            if (pred.has_value())
                                set_marking_mode(mark_id, mark, col_ids[pred.value()], c_id);
                            if (mark) 
                                set_id(col_span_start, mark_id);
                            
                            for (size_t l = col_span_start + 1; l < col_span_end; ++l) {
                                // Hit existing run
                                if (col_run_bv[l]) {
                                    set_marking_mode(mark_id, mark, col_ids[l], c_id);
                                    if (mark) 
                                        set_id(l, mark_id);
                                }
                            }

                            if (col_span_end < n && !col_run_bv[col_span_end]) {
                                set_id(col_span_end, 0);
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

    size_t save(string filename) {
        size_t written_bytes = 0;

        written_bytes += serialize_col_runs(col_runs, filename);

        std::ofstream out_ids(filename + ".col_ids");
        written_bytes += write_vec(col_run_ids, out_ids);
        return written_bytes;
    }

    /* Space saving option to mod the values before saving */
    size_t save(string filename, int id_bits) {
        assert (id_bits <= sizeof(ulint) * 8);
        ulint id_max = bit_max(id_bits);
        ulint id_bytes = bits_to_bytes(id_bits); // bytes needed to store id

        size_t written_bytes = 0;

        written_bytes += serialize_col_runs(col_runs, filename);

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

private:
    typedef struct {
        ulint interval;
        ulint offset;
        ulint height;
    } range;

    struct bitvec {
        bit_vector bv;
        size_t set_count;
        size_t max_search;

        bitvec() : bv(), set_count(0), max_search(0) {}

        bitvec(size_t size, size_t search_bound) : bv(size), set_count(0), max_search(search_bound) {}

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

        std::optional<size_t> predecessor(size_t idx) const {
            for (size_t i = idx + 1; i > 0 && i >= idx - max_search + 1; --i) {
                if (bv[i - 1]) {
                    return i - 1;
                }
            }
            return std::nullopt;
        }
    };

    FL_t& tbl;
    ulint n = 0;
    ulint r = 0;
    Mode mode = Mode::Default;
    Overlap overlap = Overlap::Append;

    bit_vector col_runs;
    vector<ulint> col_run_ids;

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

    inline void set_marking_mode(ulint& mark_id, bool& mark, ulint existing_id, ulint col_id) {
        if (existing_id > 0) {
            switch (overlap) {
                case Overlap::Append:
                    mark = false;
                    break;
                case Overlap::Split:
                    mark = true;
                    mark_id = 1;
                    break;
            }
        } else {
            mark = true;
            mark_id = col_id;
        }
    }

    void add_col_run_id(ulint id) {
        col_run_ids.push_back((id > 1) ? id - 1 : 0);
    }

    // Finds col runs and computes col_run_ids and col_runs
    void find_col_runs(std::unordered_map<ulint, ulint> &col_ids, bitvec &col_run_bv, int split_rate) {
        #ifdef PRINT_STATS
        size_t col_chars = 0;
        size_t col_bwt_runs = 0;
        size_t num_col_bits = 0;
        size_t last_idx = 0;
        #endif

        if (col_run_bv.size() == 0) {
            return;
        }

        col_runs = bit_vector(n, 0);
        bit_select col_run_select(&col_run_bv.get_bv());

        sdsl::sd_vector<> bwt_runs = tbl.get_L_heads();
        sd_select bwt_run_select(&bwt_runs);

        size_t curr_bwt_pos = 0;
        size_t col_iter = 1;
        size_t next_col_pos = col_run_select(col_iter);  
        size_t curr_id = 0;
        size_t prev_id = 0;

        auto process_col_run = [&]() {
            if (curr_id > 1 || prev_id > 1) {
                col_runs[next_col_pos] = 1;
                col_run_ids.push_back(curr_id);

                #ifdef PRINT_STATS
                ++col_bwt_runs;
                num_col_bits += curr_id > 1;

                if (col_iter > 1 && col_ids[last_idx] > 1) {
                    col_chars += curr_bwt_pos - last_idx;
                }
                last_idx = curr_bwt_pos;
                #endif
            }
        };

        for (size_t i = 1; i <= tbl.runs(); ++i) {
            curr_bwt_pos = bwt_run_select(i);

            ulint curr_id = 0;
            while (curr_bwt_pos >= next_col_pos) {
                curr_id = col_run_ids[col_iter - 1];
                if (next_col_pos < curr_bwt_pos) {
                    process_col_run();
                }

                ++col_iter;
                next_col_pos = (col_iter <= col_run_bv.set_bits()) ? col_run_select(col_iter) : tbl.size();
            }
            col_runs[curr_bwt_pos] = 1;
            col_run_ids.push_back(curr_id);
            prev_id = curr_id;

            #ifdef PRINT_STATS
            ++col_bwt_runs;
            num_col_bits += curr_id > 1;
            last_idx = curr_bwt_pos;
            #endif
        }
        while(col_iter <= col_run_bv.set_bits()) {
            curr_id = col_run_ids[col_iter - 1];
            process_col_run();
            ++col_iter;
            next_col_pos = (col_iter <= col_run_bv.set_bits()) ? col_run_select(col_iter) : tbl.size();
        }

        #ifdef PRINT_STATS
        cout << "Col ids: " << num_col_bits << std::endl;
        cout << "Col runs: " << col_bwt_runs << std::endl;
        cout << "Col chars: " << col_chars << std::endl;
        #endif
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