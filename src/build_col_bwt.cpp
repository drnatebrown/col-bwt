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

#include <common.hpp>
#include <col_bwt.hpp>

#include <sdsl/io.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <malloc_count.h>

int main(int argc, char *const argv[])
{
    Args args;
    parseArgs(argc, argv, args);

    std::cout << "Building Col BWT supporting Co-lineary Statistics: " << std::endl;
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    std::string col_ids_fname = args.filename + ".col_ids";
    std::ifstream ifs_col_ids(col_ids_fname);
    std::string col_runs_fname = args.filename + ".col_runs";
    std::ifstream ifs_col_runs(col_runs_fname);

    ifs_col_runs.seekg(0);
    sdsl::sd_vector<> col_runs;
    col_runs.load(ifs_col_runs);

    std::string bwt_fname = args.filename + ".bwt";

    std::string bwt_heads_fname = bwt_fname + ".heads";
    std::ifstream ifs_heads(bwt_heads_fname);
    std::string bwt_len_fname = bwt_fname + ".len";
    std::ifstream ifs_len(bwt_len_fname);

    ifs_heads.seekg(0);
    ifs_len.seekg(0);
    ifs_col_ids.seekg(0);
    col_bwt col_bwt(ifs_heads, ifs_len, ifs_col_ids, col_runs);

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    std::cout << "Construction Complete" << std::endl;
    std::cout << "Elapsed time (s): " << std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count() << std::endl;

    col_bwt.bwt_stats();
    col_bwt.mem_stats();

    std::cout << "Serializing" << std::endl;
    std::chrono::high_resolution_clock::time_point t_insert_mid = std::chrono::high_resolution_clock::now();

    std::string col_bwt_outfile = args.filename + col_bwt.get_file_extension();
    std::ofstream out_fp(col_bwt_outfile);
    col_bwt.serialize(out_fp);
    out_fp.close();

    t_insert_end = std::chrono::high_resolution_clock::now();

    std::cout << "Serializing Complete" << std::endl;
    std::cout << "Elapsed time (s): " << std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_mid).count() << std::endl;

    std::cout << "Done" << std::endl;
    std::cout << "Memory peak: " << malloc_count_peak() << std::endl;
    std::cout << "Total Elapsed time (s): " << std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count() << std::endl;

    return 0;
}