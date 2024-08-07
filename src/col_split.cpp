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

    col_split<> split_ds(tbl);
    split_ds.split(match_lens, match_pos, args.N);

    t_insert_end = std::chrono::high_resolution_clock::now();
    std::cout << "Splitting Complete" << std::endl;
    std::cout << "Elapsed time (s): " << std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_mid).count() << std::endl;

    std::cout << "Serializing COL runs bitvector and IDs" << std::endl;
    t_insert_mid = std::chrono::high_resolution_clock::now();

    split_ds.save(args.filename);

    t_insert_end = std::chrono::high_resolution_clock::now();
    std::cout << "Serializing Complete" << std::endl;
    std::cout << "Elapsed time (s): " << std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_mid).count() << std::endl;

    std::cout << "Done" << std::endl;
    std::cout << "Memory peak: " << malloc_count_peak() << std::endl;
    std::cout << "Total Elapsed time (s): " << std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count() << std::endl;

    return 0;
}