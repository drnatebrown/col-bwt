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
#include <FL_table.hpp>
#include <r_index.hpp>

#include <sdsl/io.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <malloc_count.h>

int main(int argc, char *const argv[])
{
    Timer timer;
    Args args;
    parseArgs(argc, argv, args);

    message("Building move structure supporting FL mapping: ");
    timer.start();

    std::string bwt_fname = args.filename + ".bwt";

    std::string bwt_heads_fname = bwt_fname + ".heads";
    std::ifstream ifs_heads(bwt_heads_fname);
    std::string bwt_len_fname = bwt_fname + ".len";
    std::ifstream ifs_len(bwt_len_fname);

    ifs_heads.seekg(0);
    ifs_len.seekg(0);
    FL_table tbl(ifs_heads, ifs_len);

    timer.end();

    tbl.bwt_stats();
    tbl.mem_stats();

    submessage("Construction Complete");
    timer.startTime();

    message("Serializing: ");
    timer.mid();

    std::string tbl_outfile = args.filename + tbl.get_file_extension();
    std::ofstream out_fp(tbl_outfile);
    tbl.serialize(out_fp);
    out_fp.close();

    timer.end();
    submessage("Serialization Complete");
    timer.midTime();

    message("Done", false);
    mem_peak();
    timer.startTime();

    return 0;
}