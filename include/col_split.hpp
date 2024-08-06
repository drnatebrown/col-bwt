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

#define PRINT_STATS 1

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
        bit_vector col_run_bv(tbl.size());
        vector<ulint> col_ids(tbl.size(), 0);
        status();

        status("FL Stepping");
        size_t run_start = 0;
        size_t c_id = 2;
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
                                col_run_bv[l] = 0; // Erase any colored runs here; they have overlaps
                            }
                            else {
                                col_ids[l] = c_id;
                            }

                            if (j % split_rate == 0 && col_ids[l] > 1 && last_id != col_ids[l]) {
                                switch (m) {
                                case Mode::Tunneled:
                                    if (FL_ranges.size() == 1) {
                                        col_run_bv[l] = 1;
                                    }
                                    break;
                                case Mode::RunAligned:
                                    if (rg.offset == 0 && tbl.get_length(rg.interval) == rg.height) {
                                        col_run_bv[l] = 1;
                                    }
                                    break;
                                default: // Mode::Default
                                    col_run_bv[l] = 1;
                                    break;
                                }
                            }

                            last_id = col_ids[l];
                        }
                        // Mark end of run if either the run was not overlapping or if the run after the overlap has a unique id
                        if (j % split_rate == 0 && col_span_end < tbl.size() && col_ids[col_span_end] > 0 && last_id != col_ids[col_span_end]) {
                            col_run_bv[col_span_end] = 1;
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

        status("Creatings IDs vector");
        #ifdef PRINT_STATS
        size_t col_chars = 0;
        size_t num_col_runs = 0;
        bool run = false;
        #endif
        for (size_t i = 0; i < tbl.size(); ++i) {
            if (col_run_bv[i] == 1) {
                col_run_ids.push_back((col_ids[i] > 1) ? col_ids[i] - 1 : 0);
            }

            #ifdef PRINT_STATS
            if (col_run_bv[i]) {
                if (col_ids[i] > 1) {
                    run = true;
                    ++num_col_runs;
                }
                else {
                    run = false;
                }
            }

            if (run) {
                ++col_chars;
            }
            #endif
        }
        status();

        status("Creating runs bitvector");
        col_runs = sd_vector(col_run_bv);
        status();

        #ifdef PRINT_STATS
        cout << " Col runs: " << num_col_runs << std::endl;
        cout << "Col chars: " << col_chars << std::endl;
        #endif
    }

    void save(string filename, int id_bits = RW_BYTES * 8) {
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
    FL_t& tbl;
    sd_vector<> col_runs;
    vector<ulint> col_run_ids;

    typedef struct {
        ulint interval;
        ulint offset;
        ulint height;
    } range;

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