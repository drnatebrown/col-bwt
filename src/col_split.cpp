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
    // ulint c_id = 1;
    assert(col_len.size() == col_pos.size());
    for (size_t i = 0; i < col_pos.size(); ++i) {
        ulint c_pos = col_pos[i];
        ulint c_len = col_len[i];
        string c_str = "";
        // add assert about the actual string
        for (size_t j = c_pos - N; j < c_pos; ++j) {
            // string new_string = "";
            // new_string += ridx.f_char(j);
            // ulint idx = ridx.FL(j);

            // assert(ridx.FL(j) == tbl.FL(ridx.get_run(j), j - ridx.run_start(ridx.get_run(j))).first);
            // ulint temp = ridx.FL(j);
            // size_t test_c = 0;
            // ulint run;
            // ulint offset;
            // for (size_t n = 0; n < tbl.runs(); ++n) {
            //     if (j < test_c + tbl.get_length(n)) {
            //         run = n;
            //         offset = j - test_c;
            //         break;
            //     }
            //     test_c += tbl.get_length(n);
            // }   
            // std::pair<ulint, ulint> jump = tbl.FL(run, offset);
            // ulint tbl_idx = tbl.get(jump.first).idx + jump.second;
            // if (temp != tbl_idx) {
            //     cerr << "ERROR: FL does not match at j=" << j << " with temp=" << temp << " and tbl_idx=" << tbl_idx << "\n";
            // }
            // else {
            //     cerr << "MATCH: FL does match at j=" << j << " with temp=" << temp << " and tbl_idx=" << tbl_idx << "\n";
            // }

            string new_string = "";
            ulint idx = j;
            for (size_t k = 0; k < c_len; ++k) {
                new_string += ridx.f_char(idx);
                idx = ridx.FL(idx);
                if (!covered_bv[idx]) ++colored_chars;
                covered_bv[idx] = 1;
                // docs[idx].push_back(c_id);
                if(std::find(docs[idx].begin(), docs[idx].end(), i) != docs[idx].end()) {
                    cerr << "DUP FOUND at l=" << idx << " with c_id=" << i << "\n";
                }
                else {
                    docs[idx].push_back(i);
                }
                assert(ridx.FL(idx) == tbl.FL(ridx.get_run(idx), idx - ridx.run_start(ridx.get_run(idx))).first);
                // idx = ridx.FL(idx);
            }
            if (c_str.length() != 0 && c_str != new_string) {
                cerr << "ERROR: Strings do not match for c_id: " << i << "\n";
                cerr << "for range: [" << c_pos << ", " << c_pos + N << "] and j=" << (j - 1) << "," << j << "\n";
                cerr << c_str << "\n";
                cerr << new_string << "\n";
                // exit(1);
            }
            c_str = new_string;
        }
        // ++c_id;
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
    cout << "Finished writing run vector\n";

    vector<vector<ulint>> docs(tbl.size());
    bv_t covered_bv(tbl.size());

    cout << "Finished initializing doc vecs\n";

    ulint colored_chars = 0;
    ulint colored_rows = 0;
    ulint total_width = 0;
    // ulint divided_rows = 0;

    ulint half_tunnels = 0;
    double max_consecutive_percent = 0;
    double avg_consecutive_percent = 0;
    ulint quarter_avg = 0;
    ulint half_avg = 0;

    ulint c_id = 1;
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

            for (size_t j = 0; j < *c_len; ++j) {
                fragments[j] = FL_ranges.size();
                vector<range> next_ranges;
                for (size_t k = 0; k < FL_ranges.size(); ++k) {
                    rg = FL_ranges[k];
                    ulint idx = tbl.get(rg.interval).idx;

                    ++colored_rows;
                    // divided_rows += (idx + rg.offset + rg.height < n && !run_bv[idx + rg.offset + rg.height]);

                    for (size_t l = idx + rg.offset; l < idx + rg.offset + rg.height; ++l) {
                        if (!covered_bv[l]) ++colored_chars;
                        covered_bv[l] = 1;
                        if(std::find(docs[l].begin(), docs[l].end(), c_id) != docs[l].end()) {
                            cerr << "DUP FOUND at l=" << l << " with c_id=" << c_id << "\n";
                        }
                        else {
                            docs[l].push_back(c_id);
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
            ulint consecutive = 0;
            for (size_t i = 0; i < *c_len; ++i) {
                if (fragments[i] == 1) {
                    ++consecutive;
                }
                else {
                    break;
                }
            }
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

    ulint curr_r_len = 0;
    ulint bwt_runs = 0;
    ulint head_runs = 0;
    ulint tail_runs = 0;
    ulint bwt_max_doc_percentage_sum = 0;
    std::unordered_map<ulint, ulint> bwt_run_frequencies;

    ulint col_runs = 0;
    ulint col_run_len = 0;
    ulint total_docs = 0;
    ulint col_max_doc_percentage_sum = 0;
    std::unordered_map<ulint, ulint> col_run_frequencies;

    ulint uniq_runs = 0;
    ulint uniq_runs_lens = 0;
    ulint uniq_count = 0;
    bool uniq = false;

    bool last = false;
    for (size_t i = 0; i < tbl.size(); ++i) {
        if (covered_bv[i]) {
            // Intersect BWT run-head
            if (run_bv[i]) {
                ++bwt_runs;
                if (docs[i].size() == 1) {
                    uniq_count += uniq_runs_lens;
                    ++uniq_runs;
                }
            }
            else if (uniq && docs[i].size() != 1) {
                uniq_count += uniq_runs_lens;
                uniq_runs_lens = 0;
            }

            // New Run
            if (!last) {
                // Corresponds to start of BWT run
                curr_r_len = 0;
                col_run_len = 0;
                uniq_runs_lens = 0;
                if (!run_bv[i]) {
                    ++head_runs;
                    if (docs[i].size() == 1) ++uniq_runs_lens;
                }
                ++col_runs;
            }
            else {
                // New BWT Run
                if (run_bv[i]) {
                    ulint max_freq = 0;
                    for (const auto& pair : bwt_run_frequencies) {
                        max_freq = std::max(max_freq, pair.second);
                    }
                    bwt_max_doc_percentage_sum += (double(max_freq) / curr_r_len)*100;
                    curr_r_len = 0;
                    uniq_runs_lens = 0;
                    bwt_run_frequencies.clear();
                }
            }

            ++col_run_len;
            ++curr_r_len;
            if (docs[i].size() == 1) ++uniq_runs_lens;
            for(const auto& doc : docs[i]) {
                bwt_run_frequencies[doc]++;
                col_run_frequencies[doc]++;
            }
        }
        // End of COL run
        else if (!covered_bv[i] && last) {
            if (!run_bv[i]) {
                ++tail_runs;
            }
            total_docs += col_run_frequencies.size();

            ulint max_freq = 0;
            for (const auto& pair : bwt_run_frequencies) {
                max_freq = std::max(max_freq, pair.second);
            }
            bwt_max_doc_percentage_sum += (double(max_freq) / curr_r_len)*100;
            curr_r_len = 0;
            bwt_run_frequencies.clear();

            max_freq = 0;
            for (const auto& pair : col_run_frequencies) {
                max_freq = std::max(max_freq, pair.second);
            }
            col_max_doc_percentage_sum += (double(max_freq) / col_run_len)*100;
            col_run_len = 0;
            col_run_frequencies.clear();
        }
        last = covered_bv[i];
        uniq = docs[i].size() == 1;
    }

    // ulint count_bv = 0;
    // ulint count_grouped = 0;
    // ulint curr_len = 0;
    // ulint curr_bwt_len = 0;
    // ulint total_bwt_len = 0;

    // ulint others_added = 0;
    // ulint extra_runs = 0;
    // ulint bwt_covered = 0;

    // ulint total_len = 0;
    // ulint total_docs = 0;
    // ulint max_doc_percentage_sum = 0;
    // ulint max_doc_percentage_sum_bwt = 0;
    // std::unordered_map<ulint, ulint> doc_frequencies;
    // std::unordered_map<ulint, ulint> doc_frequencies_bwt;

    // bool last = false;
    // for(size_t i = 0; i < covered_bv.size(); ++i) {
    //     count_bv += covered_bv[i];
    //     if (covered_bv[i]) {
    //         if (!last) {
    //             ++count_grouped;
    //             if (!run_bv[i]) {
    //                 ++extra_runs;
    //             }
    //         }
    //         else if (!last && run_bv[i]) {
    //             curr_bwt_len = 0;
    //         }
    //         else if (last && run_bv[i]) {
    //             total_bwt_len += curr_bwt_len;
    //             bwt_covered++;

    //             total_docs += doc_frequencies.size();

    //             ulint max_freq = 0;
    //             for (const auto& pair : doc_frequencies) {
    //                 max_freq = std::max(max_freq, pair.second);
    //             }
    //             max_doc_percentage_sum += (double(max_freq) / curr_len)*100;

    //             curr_len = 0;
    //             doc_frequencies.clear();
    //         }
    //         ++curr_len;
    //         for(const auto& doc : docs[i]) {
    //             doc_frequencies[doc]++;
    //             doc_frequencies_bwt[doc]++;
    //         }
    //     }
    //     else if (!covered_bv[i] && last) {
    //         others_added += !run_bv[i];
    //         total_len += curr_len;
    //         total_docs += doc_frequencies.size();

    //         ulint max_freq = 0;
    //         for (const auto& pair : doc_frequencies) {
    //             max_freq = std::max(max_freq, pair.second);
    //         }
    //         max_doc_percentage_sum += (double(max_freq) / curr_len)*100;

    //         curr_len = 0;
    //         doc_frequencies.clear();
    //     }
    //     last = covered_bv[i];
    // }
    // assert(count_bv == colored_chars);

    cout << "-----------GENERAL------------" << endl;
    cout << "            # MUMs: " << col_len.size() << endl;
    cout << "     avg MUM width: " << (double(total_width) / col_len.size()) << endl;
    cout << "                 n: " << tbl.size() << endl;
    cout << "     chars colored: " << colored_chars << endl;
    cout << "   % chars colored: " << (double(colored_chars) / (tbl.size()))*100 << endl;
    cout << "                 r: " << tbl.runs() << endl;
    cout << "                 N: " << N << endl;
    cout << "  avg MUM/COL-runs: " << (double(head_runs + bwt_runs) / col_runs) << endl;

    cout << "-----------COL RUNS-----------" << endl;
    cout << "        # COL-runs: " << col_runs << endl;
    cout << "avg COL-run height: " << (double(colored_chars) / col_runs) << endl;
    cout << "  avg BWT/COL-runs: " << (double(head_runs + bwt_runs + tail_runs) / col_runs) << endl;
    cout << "    avg docs (COL): " << (double(total_docs) / col_runs) << endl;
    cout << "    majority doc %: " << (double(col_max_doc_percentage_sum) / col_runs) << endl;

    cout << "-----------MUM RUNS-----------" << endl;
    cout << "        # MUM-runs: " << head_runs + bwt_runs << endl;
    // cout << "       # head-runs: " << head_runs << endl;
    cout << "       # tail-runs: " << tail_runs << endl;
    cout << "    total BWT-runs: " << head_runs + tail_runs + tbl.runs() << endl;
    cout << "      % r increase: " << double(head_runs + tail_runs + tbl.runs())/tbl.runs() << endl;
    cout << "        % MUM-runs: " << double(head_runs + bwt_runs)/(head_runs + tail_runs + tbl.runs())*100 << endl;
    cout << "avg MUM-run height: " << (double(colored_chars) / (head_runs + bwt_runs)) << endl;
    cout << "    avg docs (COL): " << (double(total_docs) / (head_runs + bwt_runs)) << endl;
    cout << "    majority doc %: " << (double(bwt_max_doc_percentage_sum) / (head_runs + bwt_runs)) << endl;

    cout << "-----------OVERLAPs-----------" << endl;
    cout << " # Nonoverlap-runs: " << uniq_runs << endl;
    cout << "total non+BWT-runs: " << uniq_runs + tail_runs + tbl.runs() << endl;
    cout << "  non-% r increase: " << double(uniq_runs + tail_runs + tbl.runs())/tbl.runs() << endl;
    cout << "        % MUM-runs: " << double(uniq_runs)/(head_runs + tail_runs + tbl.runs())*100 << endl;
    cout << "  avg uniq-run hgt: " << (double(uniq_count) / (uniq_runs)) << endl;
    cout << "     overlap chars: " << colored_chars - uniq_count << endl;
    cout << "        uniq chars: " << uniq_count << endl;
    cout << "      uniq chars %: " << (double(uniq_count)/tbl.size())*100 << endl;
    cout << " avg UNIQ/MUM-runs: " << (double(uniq_runs) / (head_runs + bwt_runs)) << endl;

    cout << "-----------TUNNELS------------" << endl;
    cout << "      half-tunnels: " << half_tunnels << endl;
    cout << "     half-tunnel %: " << (double(half_tunnels) / col_len.size())*100 << endl;
    cout << "   avg 1/2 tunnels: " << (double(half_avg) / col_len.size()) << endl;
    cout << "   avg 1/4 tunnels: " << (double(quarter_avg) / col_len.size()) << endl;
    cout << " avg consecutive %: " << (double(avg_consecutive_percent) / col_len.size())*100 << endl;
    cout << " max consecutive %: " << max_consecutive_percent*100 << endl;
    
    // cout << "            # MUMs: " << col_len.size() << endl;
    // cout << " Runs added (MUMs): " << count_grouped << endl;
    // cout << "Runs added (other): " << others_added << endl;
    // cout << "  Runs added (BWT): " << extra_runs + others_added << endl;
    // cout << "Runs added (other): " << others_added << endl;
    // cout << "      Colored runs: "  << (extra_runs + bwt_covered + others_added) << endl;
    // cout << " % bwt run colored: " << (double(extra_runs + bwt_covered + others_added) / (tbl.runs() + extra_runs + others_added))*100 << endl;
    // cout << "    Possible added: " << colored_rows << endl;
    // cout << "        Total rows: " << (tbl.runs() + count_grouped + others_added) << endl;
    // cout << "       Rows before: " << tbl.runs() << endl;
    // cout << "    % runs colored: " << (double(count_grouped) / (tbl.runs() + count_grouped + others_added))*100 << endl;
    // cout << "   % chars colored: " << (double(colored_chars) / (tbl.size()))*100 << endl;
    // cout << "avg colored height: " << (double(total_len) / count_grouped) << endl;
    // cout << "     avg MUM width: " << (double(total_width) / col_len.size()) << endl;
    // cout << "    avg color docs: " << (double(total_docs) / count_grouped) << endl;
    // cout << "avg color docs (r): " << (double(total_docs) / (extra_runs + bwt_covered)) << endl;
    // cout << "    majority doc %: " << (double(max_doc_percentage_sum) / count_grouped) << endl;

    // return run_bv;
}

// bv_t split(FL_table tbl, r_index<> ridx, vector<ulint> col_len, vector<ulint> col_pos, ulint N) {
//     assert(col_len.size() == col_pos.size());

//     vector<vector<ulint>> docs(tbl.size());

//     bv_t run_bv(tbl.size());
//     bv_t covered_bv(tbl.size());

//     // ulint added_rows = 0;
//     ulint colored_rows = 0;
//     ulint colored_chars = 0;
//     // ulint possible_rows = 0;

//     size_t cursor = 0;
//     ulint c_id = 1;
//     vector<ulint>::iterator c_pos = col_pos.begin();
//     vector<ulint>::iterator c_len = col_len.begin();
//     for (size_t i = 0; i < tbl.runs(); ++i) {
//         run_bv[cursor] = 1;
//         ulint run_start = cursor;
//         ulint run_len = tbl.get_length(i);

//         while (c_pos != col_pos.end() && *c_pos < run_start + run_len) {
//             // cout << i << " " << *c_pos << " " << run_start << " " << *c_len << endl;
//             // possible_rows += *c_len;
//             // ulint interval = i;
//             // ulint offset = *c_pos - run_start;
//             // ulint height = N;

//             // vector<std::pair<ulint, ulint>> col_splits;
//             // vector<ulint> split_lens;

//             // // The first FL steps reach the start of the COLs; they may diverge to different runs
//             // while (height > 0) {
//             //     // auto [inter, off] = tbl.FL(interval, offset);
//             //     // if (col_splits.size() > 0 &&inter )
//             //     col_splits.push_back(tbl.FL(interval, offset));
//             //     assert(tbl.get(col_splits.back().first).idx + col_splits.back().second == ridx.FL(*c_pos + N - height));

//             //     if (offset + height > tbl.get_length(interval))  {
//             //         ulint covered = tbl.get_length(interval) - offset;
//             //         split_lens.push_back(covered);
//             //         height -= covered;
//             //         offset = 0;

//             //     }
//             //     else {
//             //         split_lens.push_back(height);
//             //         height = 0;
//             //     }
//             //     ++interval;
//             // }

//             range rg = {i, *c_pos - run_start, N};
//             vector<range> FL_ranges = FL_range(rg, tbl, ridx);

//             for (size_t j = 0; j < *c_len; ++j) {
//                 vector<range> next_ranges;
//                 for (size_t k = 0; k < FL_ranges.size(); ++k) {
//                     rg = FL_ranges[k];
//                     ulint idx = tbl.get(rg.interval).idx;

//                     for (size_t l = idx + rg.offset; l < idx + rg.offset + rg.height; ++l) {
//                         if (!covered_bv[l]) ++colored_chars;
//                         covered_bv[l] = 1;
//                         if(std::find(docs[l].begin(), docs[l].end(), c_id) != docs[l].end()) {
//                             cerr << "DUP FOUND at l=" << l << " with c_id=" << c_id << "\n";
//                         }
//                         docs[l].push_back(c_id);
//                     }
//                     vector<range> curr_ranges = FL_range(rg, tbl, ridx);
//                     next_ranges.insert(next_ranges.end(), curr_ranges.begin(), curr_ranges.end());
//                 }
//                 FL_ranges = next_ranges;
//             }
//             301579
//             310084
//             // OLD BETTER
//             // std::pair<ulint, ulint> next = tbl.FL(interval, offset);
//             // interval = next.first;
//             // offset = next.second;

//             // ulint ri_idx = ridx.FL(*c_pos);
//             // ulint ri_idx_end = ridx.FL(*c_pos + N - 1);

//             // ulint ri_run = ridx.get_run(ri_idx);
//             // ulint ri_run_end = ridx.get_run(ri_idx_end);

//             // for(size_t j = 0; j < col_splits.size(); ++j) {
//             //     interval = col_splits[j].first;
//             //     offset = col_splits[j].second;
//             //     height = split_lens[j];
//             //     for (size_t k = 0; k < *c_len; ++k) {
//             //         ulint idx = tbl.get(interval).idx;

//             //         // added_rows += (!run_bv[idx + offset]) + (!run_bv[idx + offset + height]);
//             //         colored_rows += 1;

//             //         // run_bv[idx + offset] = 1;
//             //         // if (idx + offset + split_lens[j] < tbl.size()) run_bv[idx + offset + height] = 1;

//             //         for (size_t l = idx + offset; l < idx + offset + height; ++l) {
//             //             if (!covered_bv[l]) ++colored_chars;
//             //             covered_bv[l] = 1;
//             //             // if(std::find(docs[l].begin(), docs[l].end(), c_id) != docs[l].end()) cerr << "DUP FOUND\n";
//             //             // cout << l << " " << c_id << endl;
//             //             if(std::find(docs[l].begin(), docs[l].end(), c_id) != docs[l].end()) {
//             //                 cerr << "DUP FOUND at l=" << l << " with c_id=" << c_id << "\n";
//             //             }
//             //             docs[l].push_back(c_id);
//             //         }

//             //         std::pair<ulint, ulint> next = tbl.FL(interval, offset);
//             //         assert(tbl.get(next.first).idx + next.second == ridx.FL(tbl.get(interval).idx + offset));
//             //         interval = next.first;
//             //         offset = next.second;
//             //     }
//             // }

//             // OLD STUPID
//             // for(size_t j = 0; j < *c_len; ++j) {
//             //     ulint idx = tbl.get(interval).idx;

//             //     added_rows += (!run_bv[idx + offset]) + (!run_bv[idx + offset + N]);
//             //     colored_rows += 1;

//             //     run_bv[idx + offset] = 1;
//             //     if (idx + offset + N < tbl.size()) run_bv[idx + offset + N] = 1;
//             //     for (size_t k = idx + offset; k < idx + offset + N; ++k) {
//             //         if (!covered_bv[k]) ++colored_chars;
//             //         covered_bv[k] = 1;
//             //     }

//             //     next = tbl.FL(interval, offset);
//             //     interval = next.first;
//             //     offset = next.second;
//             // }
//             ++c_pos;
//             ++c_len;
//             ++c_id;
//         }
//         cursor = run_start + run_len;
//         printProgressBar(i + 1, tbl.runs());
//     }

//     ulint count_bv = 0;
//     ulint count_grouped = 0;
//     ulint curr_len = 0;
//     ulint others_added = 0;

//     ulint total_len = 0;
//     ulint total_docs = 0;
//     std::unordered_map<ulint, ulint> doc_frequencies;

//     std::ofstream fs_docs("temp.col_docs");

//     bool last = false;
//     for(size_t i = 0; i < covered_bv.size(); ++i) {
//         count_bv += covered_bv[i];
//         if (covered_bv[i]) {
//             if (!last) ++count_grouped;
//             ++curr_len;
//             // curr_docs.insert(docs[i]);
//             for(const auto& doc : docs[i]) {
//                 doc_frequencies[doc]++;
//             }
//         }
//         else if (!covered_bv[i] && last) {
//             others_added += !run_bv[i];
//             total_len += curr_len;
//             total_docs += doc_frequencies.size();

//             // fs_docs << curr_len << "," << curr_docs.size() << endl;
//             fs_docs << curr_len << "\t" << doc_frequencies.size() << "\t";
//             for (const auto& pair : doc_frequencies) {
//                 fs_docs << pair.second << ",";
//             }
//             fs_docs << endl;

//             curr_len = 0;
//             // curr_docs.clear();
//             fs_docs.clear();
//         }
//         last = covered_bv[i];
//     }
//     assert(count_bv == colored_chars);

//     // cout << "        Rows added: " << added_rows << endl;
//     // cout << "    Possible added: " << possible_rows*2 << endl;
//     // cout << "      Rows colored: " << colored_rows << endl;
//     cout << "       num of COLS: " << col_len.size() << endl;
//     cout << "  Rows added (col): " << count_grouped << endl;
//     cout << "  Rows added (oth): " << others_added << endl;
//     cout << "        Total rows: " << (tbl.runs() + count_grouped + others_added) << endl;
//     cout << "       Rows before: " << tbl.runs() << endl;
//     cout << "    % runs colored: " << (double(count_grouped) / (tbl.runs() + count_grouped + others_added))*100 << endl;
//     cout << "   % chars colored: " << (double(colored_chars) / (tbl.size()))*100 << endl;
//     cout << "avg colored height: " << (double(total_len) / count_grouped) << endl;
//     cout << "      avg num docs: " << (double(total_docs) / count_grouped) << endl;
//     // cout << " avg colored width: " << (double(colored_rows) / col_len.size()) << endl;

//     return run_bv;
// }

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