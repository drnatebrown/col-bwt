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

template <class FL_t = FL_table>
class col_split {
public:
    enum class Mode {
        Default,
        Tunneled,
        RunAligned
    };

    col_split() {}

    col_split(FL_t& t) : tbl(t) {}

    void split(vector<ulint> col_len, vector<ulint> col_pos, ulint N, Mode m = Mode::Default, int split_rate = 1) {
        assert(col_len.size() == col_pos.size());

        status("Initializing vectors");
        bitvec col_run_bv(tbl.size());
        vector<ulint> col_ids(tbl.size(), 0);
        status();

        status("FL Stepping");
        size_t run_start = 0;
        size_t c_id = 2; // 0 denotes to id, 1 denotes overlap
        vector<ulint>::iterator c_pos = col_pos.begin();
        vector<ulint>::iterator c_len = col_len.begin();
        for (size_t i = 0; i < tbl.runs(); ++i) {
            ulint run_len = tbl.get_length(i);

            while (c_pos != col_pos.end() && *c_pos - N < run_start + run_len) {
                range rg = {i, *c_pos - N - run_start, N};
                vector<range> FL_ranges = FL_range(rg);

                for (size_t j = 0; j < *c_len; ++j) {
                    vector<range> next_ranges;
                    for (size_t k = 0; k < FL_ranges.size(); ++k) {
                        rg = FL_ranges[k];
                        ulint idx = tbl.get_idx(rg.interval);
                        ulint col_span_start = idx + rg.offset;
                        ulint col_span_end = idx + rg.offset + rg.height;

                        ulint last_id = (col_span_start > 0) ? col_ids[col_span_start - 1] : 0;
                        for (size_t l = col_span_start; l < col_span_end; ++l) {
                            if (col_ids[l] >= 1) {
                                col_ids[l] = 1;
                                col_run_bv.unset(l); // Erase any colored runs here; they have overlaps
                            }
                            else {
                                col_ids[l] = c_id;
                            }

                            // Split if split rate is met and new run or we are splitting another run due to overlap
                            if (last_id != col_ids[l] && ((j % split_rate == 0 && col_ids[l] == c_id) || col_ids[l] == 1 || last_id == 1)) {
                                switch (m) {
                                case Mode::Tunneled:
                                    if (FL_ranges.size() == 1) {
                                        col_run_bv.set(l);
                                    }
                                    break;
                                case Mode::RunAligned:
                                    if (rg.offset == 0 && tbl.get_length(rg.interval) == rg.height) {
                                        col_run_bv.set(l);
                                    }
                                    break;
                                default: // Mode::Default
                                    col_run_bv.set(l);
                                    break;
                                }
                            }
                            else if (last_id == col_ids[l]) {
                                col_run_bv.unset(l);
                            }

                            last_id = col_ids[l];
                        }
                        
                        if (col_span_end < tbl.size()) {
                            // End the current run if:
                            // i) we are at the end of a run of overlapping characters and the next is not overlapping
                            // ii) we are at the end of a run of non-overlapping characters and the next id does not match
                            if ((last_id > 1 || col_ids[col_span_end] > 0) && last_id != col_ids[col_span_end]) {
                                col_run_bv.set(col_span_end);
                            }
                            // If we can continue the run, do so
                            else if (last_id == col_ids[col_span_end]) {
                                col_run_bv.unset(col_span_end);
                            }
                        }

                        // TODO: Paint non-overlapping run aligned even if split rate not met

                        vector<range> curr_ranges = FL_range(rg);
                        next_ranges.insert(next_ranges.end(), curr_ranges.begin(), curr_ranges.end());
                    }
                    FL_ranges = next_ranges;
                }
                ++c_pos;
                ++c_len;
                ++c_id;
            }
            run_start += run_len;
        }
        status();

        status("Creating runs bitvector");
        col_runs = sd_vector(col_run_bv.get_bv());
        sd_select col_select = sd_select(&col_runs);
        status();

        status("Creatings IDs vector");
        #ifdef PRINT_STATS
        size_t col_chars = 0;
        size_t num_col_runs = 0;
        bool run = false;
        size_t last_idx = 0;
        #endif
        for (size_t i = 1; i <= col_run_bv.set_bits(); ++i) {
            size_t idx = col_select(i);
            col_run_ids.push_back((col_ids[idx] > 1) ? col_ids[idx] - 1 : 0);

            #ifdef PRINT_STATS
            if (run) {
                col_chars += idx - last_idx;
            }

            if (col_ids[idx] > 1) {
                run = true;
                last_idx = idx;
                ++num_col_runs;
            }
            else {
                run = false;
            }
            #endif
        }
        status();

        #ifdef PRINT_STATS
        cout << "Col runs: " << num_col_runs << std::endl;
        cout << "Col chars: " << col_chars << std::endl;
        #endif
    }

    // void save(string filename) {
    //     std::ofstream out_runs(filename + ".col_runs");
    //     col_runs.serialize(out_runs);
    //     out_runs.close();

    //     std::ofstream out_ids(filename + ".col_ids");
    //     for (size_t i = 0; i < col_run_ids.size(); ++i) {
    //         out_ids.write(reinterpret_cast<const char*>(&col_run_ids[i]), sizeof(RW_BYTES));
    //     }
    // }

    /* Space saving option to mod the values before saving */
    void save(string filename, int id_bits=ID_BITS) {
        assert (id_bits <= sizeof(ulint) * 8);
        ulint id_max = 1 << id_bits;
        ulint id_bytes = (id_bits + 7) / 8; // bytes needed to store id

        std::ofstream out_runs(filename + ".col_runs");
        col_runs.serialize(out_runs);
        out_runs.close();

        std::ofstream out_ids(filename + ".col_ids");
        for (size_t i = 0; i < col_run_ids.size(); ++i) {
            size_t id = col_run_ids[i] % id_max;
            out_ids.write(reinterpret_cast<const char*>(&id), id_bytes);
        }
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

    FL_t& tbl;
    sd_vector<> col_runs;
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
};

#endif /* end of include guard: _COL_SPLIT_HH */