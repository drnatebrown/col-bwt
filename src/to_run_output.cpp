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
   \file cil_split.cpp
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
#include <col_split.hpp>

#include <iostream>
#include <sdsl/io.hpp>
#include <sdsl/int_vector.hpp>
#include <malloc_count.h>
#include <sys/stat.h>

// int main(int argc, char *const argv[])
// {
//     Timer timer;
//     Args args;
//     parseArgs(argc, argv, args);

//     message("Creating Movi useable output: ");
//     timer.start();

//     std::string col_ids_fname = args.filename + ".col_ids";
//     std::ifstream ifs_col_ids(col_ids_fname);
//     std::string col_runs_fname = args.filename + ".col_runs";
//     std::ifstream ifs_col_runs(col_runs_fname);

//     status("Loading col runs");
//     sdsl::sd_vector<> col_runs;
//     col_runs.load(ifs_col_runs);
//     status();

//     status("Loading FL table");
//     FL_table tbl;
//     std::string filename_tbl = args.filename + tbl.get_file_extension();
//     ifstream fs_tbl(filename_tbl);
//     tbl.load(fs_tbl);
//     sdsl::sd_vector<> bwt_runs = tbl.get_L_heads();
//     status();

//     status("Loading col ids");
//     vector<ulint> col_run_ids;
//     ifs_col_ids.clear();
//     ifs_col_ids.seekg(0);

//     size_t curr_id = 0;
//     while (ifs_col_ids.read((char *)&curr_id, RW_BYTES))
//     {
//         col_run_ids.push_back(curr_id);
//     }
//     status();

//     status("Creating new col runs/ids");
//     bit_vector new_bv(tbl.size(), 0);
//     vector<ulint> new_ids;

//     size_t curr_col_run = 0;
//     curr_id = 0;
//     for (size_t i = 0; i < tbl.size(); ++i) {
//         if (col_runs[i]) {
//             curr_id = col_run_ids[curr_col_run++];
//         }

//         if (col_runs[i] || bwt_runs[i]) {
//             new_bv[i] = 1;
//             new_ids.push_back(curr_id);
//         }
//     }
//     status();

//     status("Saving new col runs/ids");
//     std::ofstream ofs_col_runs(args.filename + ".movi_col_runs");
//     std::ofstream ofs_col_ids(args.filename + ".movi_col_ids");

//     size_t max = bit_max(8);
//     new_bv.serialize(ofs_col_runs);
//     for (size_t i = 0; i < new_ids.size(); ++i) {
//         size_t c_id = new_ids[i];
//         if (c_id >= max) {
//             c_id = (c_id % (max - 1)) + 1;
//         }
//         ofs_col_ids.write(reinterpret_cast<const char*>(&c_id), 1);
//     }
//     status();

//     return 0;
// }

int main(int argc, char *const argv[]) {
    Args args;
    parseArgs(argc, argv, args);

    if (args.long_pattern) {
        std::cout << ID_BYTES << std::endl;
        std::cout << ID_BITS << std::endl;

        Timer timer;
        message("Creating Movi useable output: ");
        timer.start();

        std::string col_ids_fname = args.filename + ".col_ids";
        std::ifstream ifs_col_ids(col_ids_fname);
        std::string col_runs_fname = args.filename + ".col_runs";
        std::ifstream ifs_col_runs(col_runs_fname);

        status("Loading col runs");
        sdsl::sd_vector<> col_runs;
        col_runs.load(ifs_col_runs);
        status();

        status("Loading FL table");
        FL_table tbl;
        std::string filename_tbl = args.filename + tbl.get_file_extension();
        ifstream fs_tbl(filename_tbl);
        tbl.load(fs_tbl);
        sdsl::sd_vector<> bwt_runs = tbl.get_L_heads();
        status();

        status("Loading col ids");
        vector<ulint> col_run_ids;
        ifs_col_ids.clear();
        ifs_col_ids.seekg(0);

        size_t curr_id = 0;
        while (ifs_col_ids.read((char *)&curr_id, RW_BYTES))
        {
            col_run_ids.push_back(curr_id);
        }
        status();

        status("Creating new col runs/ids");
        bit_vector new_bv(tbl.size(), 0);
        vector<ulint> new_ids;

        sd_select bwt_run_select(&bwt_runs);
        sd_select col_run_select(&col_runs);

        size_t curr_bwt_pos = 0;

        size_t next_col_run = 1;
        size_t next_col_pos = col_run_select(next_col_run);
        curr_id = 0;

        for (size_t i = 1; i <= tbl.runs(); ++i) {
            curr_bwt_pos = bwt_run_select(i);
            while (curr_bwt_pos >= next_col_pos) {
                curr_id = col_run_ids[next_col_run - 1];
                if (next_col_pos < curr_bwt_pos) {
                    new_bv[next_col_pos] = 1;
                    new_ids.push_back(curr_id);
                }

                ++next_col_run;
                next_col_pos = (next_col_run <= col_run_ids.size()) ? col_run_select(next_col_run) : tbl.size();
            }
            new_bv[curr_bwt_pos] = 1;
            new_ids.push_back(curr_id);
        }
        while(next_col_run <= col_run_ids.size()) {
            new_bv[col_run_select(next_col_run)] = 1;
            new_ids.push_back(col_run_ids[next_col_run - 1]);
            ++next_col_run;
        }

        // for (size_t i = 0; i < tbl.size(); ++i) {
        //     if (col_runs[i]) {
        //         curr_id = col_run_ids[curr_col_run++];
        //     }

        //     if (col_runs[i] || bwt_runs[i]) {
        //         new_bv[i] = 1;
        //         new_ids.push_back(curr_id);
        //     }
        // }
        status();

        status("Saving new col runs/ids");
        std::ofstream ofs_col_runs(args.filename + ".movi_col_runs");
        std::ofstream ofs_col_ids(args.filename + ".movi_col_ids");

        size_t max = bit_max(8);
        new_bv.serialize(ofs_col_runs);
        for (size_t i = 0; i < new_ids.size(); ++i) {
            size_t c_id = col_id_bin(new_ids[i], max);
            ofs_col_ids.write(reinterpret_cast<const char*>(&c_id), 1);
        }
        status();
        return 0;
    }
    else {
        message("Creating BWT from RLBWT");
        std::string head_filename = args.filename + ".bwt.heads";
        std::string len_filename = args.filename + ".bwt.len";
        std::string output_filename = args.filename + ".bwt";

        std::ifstream head_file(head_filename, std::ios::binary);
        std::ifstream len_file(len_filename, std::ios::binary);
        std::ofstream output_file(output_filename, std::ios::binary);

        if (!head_file.is_open() || !len_file.is_open() || !output_file.is_open()) {
            std::cerr << "Error: Could not open one or more files." << std::endl;
            std::cerr << "head_file: " << head_filename << std::endl;
            std::cerr << "len_file: " << len_filename << std::endl;
            std::cerr << "output_file: " << output_filename << std::endl;
            return 1;
        }

        status("Creating BWT");

        char bwt_char;
        uint64_t bwt_len;
        while (head_file.get(bwt_char) && len_file.read(reinterpret_cast<char*>(&bwt_len), RW_BYTES)) {
            std::string output_string(bwt_len, bwt_char);
            output_file.write(output_string.c_str(), output_string.size());
        }
        status();

        head_file.close();
        len_file.close();
        output_file.close();
        return 0;
    }
}