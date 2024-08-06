/* build_FL - Table supporting first-to-last mapping of BWT
    Copyright (C) 2024 Nathaniel Brown
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file build_FL.cpp
   \brief Builds table supporting first-to-last mapping of BWT
   \author Nathaniel Brown
   \date 06/13/2024
*/

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

typedef sdsl::bit_vector bv_t;

typedef struct {
    ulint interval;
    ulint offset;
    ulint height;
} range;

vector<range> FL_range(range rg, FL_table& tbl, r_index<>& ridx) {
    vector<range> FL_dest_ranges;

    // FL steps may diverge to different positions even within the same run
    while (rg.height > 0) {
        const auto& [FL_interval, FL_offset] = tbl.FL(rg.interval, rg.offset);
        assert(tbl.get(FL_interval).idx + FL_offset == ridx.FL(tbl.get(rg.interval).idx + rg.offset));

        if (rg.offset + rg.height > tbl.get_length(rg.interval))  {
            ulint covered = tbl.get_length(rg.interval) - rg.offset;
            // assert(covered == ridx.run_len(ridx.get_run(tbl.get(rg.interval).idx + rg.offset)) - rg.offset);
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

bv_t split_r(FL_table tbl, r_index<> ridx, vector<ulint> col_len, vector<ulint> col_pos, ulint N) {
    bv_t run_bv(ridx.size());
    size_t cursor = 0;
    for (size_t i = 0; i < ridx.runs(); ++i) {
        run_bv[cursor] = 1;
        cursor += ridx.run_len(i);
    }
    cout << "Finished writing run vector\n";

    vector<vector<ulint>> docs(ridx.size());
    bv_t covered_bv(ridx.size());

    cout << "Finished initializing doc vecs\n";

    ulint colored_chars = 0;
    assert(col_len.size() == col_pos.size());
    for (size_t i = 0; i < col_pos.size(); ++i) {
        ulint c_pos = col_pos[i];
        ulint c_len = col_len[i];
        string c_str = "";

        for (size_t j = c_pos - N; j < c_pos; ++j) {
            string new_string = "";
            ulint idx = j;
            for (size_t k = 0; k < c_len; ++k) {
                new_string += ridx.f_char(idx);
                idx = ridx.FL(idx);
                if (!covered_bv[idx]) ++colored_chars;
                covered_bv[idx] = 1;
                assert(std::find(docs[idx].begin(), docs[idx].end(), i) != docs[idx].end());
                assert(ridx.FL(idx) == tbl.FL(ridx.get_run(idx), idx - ridx.run_start(ridx.get_run(idx))).first);
                docs[idx].push_back(i + 1);
            }
        }
        if (i % (col_pos.size()/100) == 0) cout << i << "/" << col_pos.size() << endl;
    }

    ulint count_bv = 0;
    ulint count_grouped = 0;
    ulint curr_len = 0;
    ulint others_added = 0;

    ulint total_len = 0;
    ulint total_docs = 0;
    std::unordered_map<ulint, ulint> doc_frequencies;

    std::ofstream fs_docs("temp.col_docs");

    bool last = false;
    for(size_t i = 0; i < covered_bv.size(); ++i) {
        count_bv += covered_bv[i];
        if (covered_bv[i]) {
            if (!last) ++count_grouped;
            ++curr_len;
            for(const auto& doc : docs[i]) {
                doc_frequencies[doc]++;
            }
        }
        else if (!covered_bv[i] && last) {
            others_added += !run_bv[i];
            total_len += curr_len;
            total_docs += doc_frequencies.size();

            fs_docs << curr_len << "\t" << doc_frequencies.size() << "\t";
            for (const auto& pair : doc_frequencies) {
                fs_docs << pair.second << ",";
            }
            fs_docs << endl;

            curr_len = 0;
            doc_frequencies.clear();
        }
        last = covered_bv[i];
        if (i % (ridx.size()/100) == 0) cout << i << "/" << ridx.size() << endl;
    }
    assert(count_bv == colored_chars);

    // cout << "        Rows added: " << added_rows << endl;
    // cout << "    Possible added: " << possible_rows*2 << endl;
    // cout << "      Rows colored: " << colored_rows << endl;
    cout << "       num of COLS: " << col_len.size() << endl;
    cout << "  Rows added (col): " << count_grouped << endl;
    cout << "  Rows added (oth): " << others_added << endl;
    cout << "        Total rows: " << (ridx.runs() + count_grouped + others_added) << endl;
    cout << "       Rows before: " << ridx.runs() << endl;
    cout << "    % runs colored: " << (double(count_grouped) / (ridx.runs() + count_grouped + others_added))*100 << endl;
    cout << "   % chars colored: " << (double(colored_chars) / (ridx.size()))*100 << endl;
    cout << "avg colored height: " << (double(total_len) / count_grouped) << endl;
    cout << "      avg num docs: " << (double(total_docs) / count_grouped) << endl;
    // cout << " avg colored width: " << (double(colored_rows) / col_len.size()) << endl;

    return run_bv;
}

void split(FL_table tbl, r_index<> ridx, vector<ulint> col_len, vector<ulint> col_pos, ulint N) {
    assert(col_len.size() == col_pos.size());
    bv_t run_bv(ridx.size());
    size_t cursor = 0;
    for (size_t i = 0; i < ridx.runs(); ++i) {
        run_bv[cursor] = 1;
        cursor += ridx.run_len(i);
    }

    bv_t covered_bv(tbl.size());
    bv_t covered_tunnel_bv(tbl.size());
    bv_t covered_non_bv(tbl.size());
    bv_t same_run_bv(tbl.size());
    vector<ulint> col_ids(tbl.size(), 0);

    cout << "Finished initializing vectors\n";

    ulint colored_chars = 0;
    ulint colored_rows = 0;
    ulint total_width = 0;

    ulint tunneled_runs = 0;

    ulint runny_cov = 0;
    ulint runny_num = 0;

    ulint half_tunnels = 0;
    double max_consecutive_percent = 0;
    double avg_consecutive_percent = 0;
    ulint quarter_avg = 0;
    ulint half_avg = 0;

    ulint c_id = 2;
    cursor = 0;
    vector<ulint>::iterator c_pos = col_pos.begin();
    vector<ulint>::iterator c_len = col_len.begin();
    for (size_t i = 0; i < tbl.runs(); ++i) {
        ulint run_start = cursor;
        ulint run_len = tbl.get_length(i);

        while (c_pos != col_pos.end() && *c_pos - N < run_start + run_len) {
            total_width += *c_len;
            range rg = {i, *c_pos - N - run_start, N};
            vector<range> FL_ranges = FL_range(rg, tbl, ridx);
            vector<ulint> fragments(*c_len, 0);
            ulint consecutive = 0;
            for (size_t j = 0; j < *c_len; ++j) {
                fragments[j] = FL_ranges.size();
                if (fragments[j] == 1) {
                    ++consecutive;
                }
                vector<range> next_ranges;
                for (size_t k = 0; k < FL_ranges.size(); ++k) {
                    rg = FL_ranges[k];
                    ulint idx = tbl.get(rg.interval).idx;

                    ++colored_rows;
                    
                    if (idx + rg.offset + rg.height < tbl.size()) {
                        if (run_bv[idx + rg.offset] && run_bv[idx + rg.offset + rg.height]) {
                            if (!same_run_bv[idx + rg.offset]) {
                                ++runny_num;
                                runny_cov += rg.height;
                            }
                            same_run_bv[idx + rg.offset] = 1;
                        }
                    }

                    for (size_t l = idx + rg.offset; l < idx + rg.offset + rg.height; ++l) {
                        if (col_ids[l] >= 1) {
                            col_ids[l] = 1;
                        }
                        else {
                            col_ids[l] = c_id;
                        }
                        if (!covered_bv[l]) ++colored_chars;
                        covered_bv[l] = 1;

                        if (FL_ranges.size() == 1) {
                            // ++tunnels_cov;
                            covered_tunnel_bv[l] = 1;
                        }

                        if (col_ids[l] > 1) {
                            covered_non_bv[l] = 1;
                        }
                        else
                        {
                            covered_non_bv[l] = 0;
                        }
                    }
                    vector<range> curr_ranges = FL_range(rg, tbl, ridx);
                    next_ranges.insert(next_ranges.end(), curr_ranges.begin(), curr_ranges.end());
                }
                FL_ranges = next_ranges;
            }

            if (fragments[*c_len/2] == 1) half_tunnels++;
            quarter_avg += fragments[*c_len/4];
            half_avg += fragments[*c_len/2];
            avg_consecutive_percent += double(consecutive)/(*c_len);
            max_consecutive_percent = std::max(double(consecutive)/(*c_len), max_consecutive_percent);

            ++c_pos;
            ++c_len;
            ++c_id;
        }
        cursor = run_start + run_len;
        // if (i % (tbl.runs()/100) == 0) cout << i << "/" << tbl.runs() << endl;
    }

    cout << "Finished FL stepping\n";

    bv_t col_run_bv(tbl.size());
    bv_t non_run_bv(tbl.size());

    ulint curr_r_len = 0;
    ulint colored_runs = 0;
    ulint tail_runs = 0;
    ulint old_runs = 0;

    ulint curr_r_len_non = 0;
    ulint colored_runs_non = 0;
    ulint tail_runs_non = 0;
    ulint old_runs_non = 0;

    ulint curr_r_len_tnl = 0;
    ulint colored_runs_tnl = 0;
    ulint tail_runs_tnl = 0;
    ulint old_runs_tnl = 0;

    ulint colored_chars_non = 0;
    ulint colored_chars_tnl = 0;
    bool last = false;
    bool last_non = false;
    bool last_tnl = false;
    ulint last_id = 0;
    for (size_t i = 0; i < tbl.size(); ++i) {
        colored_chars_non += covered_non_bv[i];
        // covered_tunnel_bv[i] = covered_tunnel_bv[i] && covered_non_bv[i];
        colored_chars_tnl += covered_tunnel_bv[i];

        if (covered_bv[i]) {
            // New Run
            if (!last || col_ids[i] != last_id) {
                // Corresponds to start of BWT run
                curr_r_len = 0;
                ++colored_runs;
                if (run_bv[i]) ++old_runs;
            }
            else {
                // New BWT Run
                if (run_bv[i]) {
                    curr_r_len = 0;
                    ++colored_runs;
                    ++old_runs;
                }
            }
            ++curr_r_len;
        }
        // End of COL run
        else if (!covered_bv[i] && last) {
            ++tail_runs;
            if (run_bv[i]) ++old_runs;
            
            curr_r_len = 0;
        }

        if (covered_non_bv[i]) {
            // New Run
            if (!last_non || col_ids[i] != last_id) {
                // Corresponds to start of BWT run
                curr_r_len_non = 0;
                ++colored_runs_non;
                if (run_bv[i]) ++old_runs_non;
            }
            else {
                // New BWT Run
                if (run_bv[i]) {
                    curr_r_len_non = 0;
                    ++colored_runs_non;
                    ++old_runs_non;
                }
            }
            ++curr_r_len_non;
        }
        // End of COL run
        else if (!covered_non_bv[i] && last_non) {
            ++tail_runs_non;
            if (run_bv[i]) ++old_runs_non;
            
            curr_r_len_non = 0;
        }

        if (covered_tunnel_bv[i]) {
            // New Run
            if (!last_tnl || col_ids[i] != last_id) {
            // if (!last_tnl) {
                // Corresponds to start of BWT run
                curr_r_len_tnl = 0;
                ++colored_runs_tnl;
                if (run_bv[i]) ++old_runs_tnl;
            }
            else {
                // New BWT Run
                if (run_bv[i]) {
                    curr_r_len_tnl = 0;
                    ++colored_runs_tnl;
                    ++old_runs_tnl;
                }
            }
            ++curr_r_len_tnl;
        }
        // End of COL run
        else if (!covered_tunnel_bv[i] && last_tnl) {
            ++tail_runs_tnl;
            if (run_bv[i]) ++old_runs_tnl;
            
            curr_r_len_tnl = 0;
        }

        last = covered_bv[i];
        last_non = covered_non_bv[i];
        last_tnl = covered_tunnel_bv[i];
        last_id = col_ids[i];
    }

    cout << "-----------GENERAL------------" << endl;
    cout << "            # MUMs: " << col_len.size() << endl;
    cout << "     avg MUM width: " << (double(total_width) / col_len.size()) << endl;
    cout << "                 n: " << tbl.size() << endl;
    cout << "                 r: " << tbl.runs() << endl;
    cout << "               n/r: " << double(tbl.size())/tbl.runs() << endl;
    cout << "                 N: " << N << endl;

    ulint total_runs = (colored_runs + tail_runs) - old_runs + tbl.runs();
    cout << "-----------MUM RUNS-----------" << endl;
    cout << "        # MUM-runs: " << colored_runs << endl;
    cout << "  # Intersect-runs: " << old_runs << endl;
    cout << "       # tail-runs: " << tail_runs << endl;
    cout << "    total BWT-runs: " << total_runs << endl;
    cout << "      % r increase: " << double(total_runs)/tbl.runs() << endl;
    cout << "        % MUM-runs: " << (double(colored_runs)/total_runs)*100 << endl;
    cout << "avg MUM-run height: " << (double(colored_chars) / (colored_runs)) << endl;
    cout << "     chars colored: " << colored_chars << endl;
    cout << "   % chars colored: " << (double(colored_chars) / (tbl.size()))*100 << endl;
    cout << "               n/r: " << double(tbl.size())/total_runs << endl;

    total_runs = (colored_runs_non + tail_runs_non) - old_runs_non + tbl.runs();
    cout << "-----------NON RUNS-----------" << endl;
    cout << "        # NON-runs: " << colored_runs_non << endl;
    cout << "  # Intersect-runs: " << old_runs_non << endl;
    cout << "       # tail-runs: " << tail_runs_non << endl;
    cout << "    total BWT-runs: " << total_runs << endl;
    cout << "      % r increase: " << double(total_runs)/tbl.runs() << endl;
    cout << "        % NON-runs: " << (double(colored_runs_non)/total_runs)*100 << endl;
    cout << "avg NON-run height: " << (double(colored_chars_non) / (colored_runs_non)) << endl;
    cout << "     chars colored: " << colored_chars_non << endl;
    cout << "   % chars colored: " << (double(colored_chars_non) / (tbl.size()))*100 << endl;
    cout << "               n/r: " << double(tbl.size())/total_runs << endl;

    cout << "-----------TUNNELS------------" << endl;
    cout << "      half-tunnels: " << half_tunnels << endl;
    cout << "     half-tunnel %: " << (double(half_tunnels) / col_len.size())*100 << endl;
    cout << "   avg 1/2 tunnels: " << (double(half_avg) / col_len.size()) << endl;
    cout << "   avg 1/4 tunnels: " << (double(quarter_avg) / col_len.size()) << endl;
    cout << " avg consecutive %: " << (double(avg_consecutive_percent) / col_len.size())*100 << endl;
    cout << " max consecutive %: " << max_consecutive_percent*100 << endl;

    cout << "------------RUNNY-------------" << endl;
    cout << "        # RUN-runs: " << runny_num << endl;
    cout << "        % RUN-runs: " << (double(runny_num)/tbl.runs())*100 << endl;
    cout << "avg RUN-run height: " << (double(runny_cov) / (runny_num)) << endl;
    cout << "     chars colored: " << runny_cov << endl;
    cout << "   % chars colored: " << (double(runny_cov) / (tbl.size()))*100 << endl;

    total_runs = (colored_runs_tnl + tail_runs_tnl) - old_runs_tnl + tbl.runs();
    cout << "-----------TUN RUNS-----------" << endl;
    cout << "        # TUN-runs: " << colored_runs_tnl << endl;
    cout << "  # Intersect-runs: " << old_runs_tnl << endl;
    cout << "       # tail-runs: " << tail_runs_tnl << endl;
    cout << "    total BWT-runs: " << total_runs << endl;
    cout << "      % r increase: " << double(total_runs)/tbl.runs() << endl;
    cout << "        % TUN-runs: " << (double(colored_runs_tnl)/total_runs)*100 << endl;
    cout << "avg TUN-run height: " << (double(colored_chars_tnl) / (colored_runs_tnl)) << endl;
    cout << "     chars colored: " << colored_chars_tnl << endl;
    cout << "   % chars colored: " << (double(colored_chars_tnl) / (tbl.size()))*100 << endl;
    cout << "               n/r: " << double(tbl.size())/total_runs << endl;
}

int main(int argc, char *const argv[])
{
    Args args;
    parseArgs(argc, argv, args);

    if (!args.N) {
        std::cout << "ERROR: Invalid number of arguments\n";
        return 1;
    }

    std::cout << "Loading move structure supporting FL mapping: " << std::endl;
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    FL_table tbl;

    std::string filename_tbl = args.filename + tbl.get_file_extension();
    ifstream fs_tbl(filename_tbl);

    tbl.load(fs_tbl);

    r_index<> ridx;
    std::string filename_ridx = args.filename + ridx.get_file_extension();
    ifstream fs_ridx(filename_ridx);

    ridx.load(fs_ridx);

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    std::cout << "Load Complete" << std::endl;
    std::cout << "Elapsed time (s): " << std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count() << std::endl;

    tbl.bwt_stats();

    std::chrono::high_resolution_clock::time_point t_insert_mid = std::chrono::high_resolution_clock::now();
    std::cout << "Loading COL Positions" << std::endl;

    std::string filename_mums = args.filename + ".mums";
    std::ifstream mum_file(filename_mums, std::ios::binary | std::ios::ate);
    size_t num_mums = mum_file.tellg()/(RW_BYTES*4);

    std::cout << "There are " << num_mums << " COL positions" << std::endl;
    mum_file.seekg(0, std::ios::beg);

    std::vector<ulint> match_lens(num_mums);
    std::vector<ulint> match_pos(num_mums);

    char discard[RW_BYTES];
    for (size_t i = 0; i < num_mums; ++i) {
        mum_file.read(reinterpret_cast<char*>(&match_lens[i]), RW_BYTES);
        mum_file.read(reinterpret_cast<char*>(&match_pos[i]), RW_BYTES);

        // Discard the next two values
        mum_file.read(discard, RW_BYTES);
        mum_file.read(discard, RW_BYTES);
    }
    t_insert_end = std::chrono::high_resolution_clock::now();
    std::cout << "Load Complete" << std::endl;
    std::cout << "Elapsed time (s): " << std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_mid).count() << std::endl;

    std::cout << "Splitting runs based on COL Positions using FL Table" << std::endl;
    t_insert_mid = std::chrono::high_resolution_clock::now();

    split(tbl, ridx, match_lens, match_pos, args.N);

    t_insert_end = std::chrono::high_resolution_clock::now();
    std::cout << "Splitting Complete" << std::endl;
    std::cout << "Elapsed time (s): " << std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_mid).count() << std::endl;

    // std::cout << "Splitting runs based on COL Positions using r-index" << std::endl;
    // t_insert_mid = std::chrono::high_resolution_clock::now();

    // bv_t split_bv_2 = split_r(tbl, ridx, match_lens, match_pos, args.N);

    // t_insert_end = std::chrono::high_resolution_clock::now();
    // std::cout << "Splitting Complete" << std::endl;
    // std::cout << "Elapsed time (s): " << std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_mid).count() << std::endl;

    // std::cout << "Serializing runs bitvector" << std::endl;

    // std::string split_bv_fname = args.filename + ".rbv";
    // std::ofstream split_bv_file(split_bv_fname);
    // split_bv.serialize(split_bv_file);

    // t_insert_end = std::chrono::high_resolution_clock::now();
    // std::cout << "Serializing Complete" << std::endl;
    // std::cout << "Elapsed time (s): " << std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_mid).count() << std::endl;

    std::cout << "Done" << std::endl;
    std::cout << "Memory peak: " << malloc_count_peak() << std::endl;
    std::cout << "Total Elapsed time (s): " << std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count() << std::endl;

  return 0;
}