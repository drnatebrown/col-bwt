/* pml_query - Compute pml using Col ids
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
   \file pml_query.cpp
   \brief Compute pml using Col ids
   \author Nathaniel Brown
   \date 06/13/2024
*/

#include <cassert>
#include <common.hpp>
#include <col_bwt.hpp>
#include <cstddef>
#include <fstream>

#include <sdsl/io.hpp>
#include <sdsl/int_vector.hpp>
#include <malloc_count.h>

int main(int argc, char *const argv[])
{
    Timer timer;
    Args args;
    parseArgs(argc, argv, args);

    message("Loading BWT table supporting LF mapping: ");
    timer.start();

    col_pml tbl;
    std::string filename_tbl = args.filename + tbl.get_file_extension();
    ifstream fs_tbl(filename_tbl);
    tbl.load(fs_tbl);

    tbl.bwt_stats();

    timer.end();
    submessage("Load Complete");
    timer.startTime();

    timer.mid();
    message("Reading pattern file: ");

    std::ifstream fs_pattern(args.pattern_filename);
    std::string pattern;
    std::string line;
    std::ostringstream sequence_stream;
    while (std::getline(fs_pattern, line))
    {
        if (line.empty() || line[0] == '>')
            continue;
        sequence_stream << line;
    }
    pattern = sequence_stream.str();

    timer.end();
    submessage("Read Complete");
    timer.midTime();

    timer.mid();
    message("Computing PMLs Query: ");

    std::string pml_filename = args.pattern_filename + ".pml";
    std::string cid_filename = args.pattern_filename + ".cid";
    std::ofstream fs_pml(pml_filename);
    std::ofstream fs_cid(cid_filename);
    tbl.query_pml(pattern, fs_pml, fs_cid);

    timer.end();
    submessage("Query Complete");
    timer.midTime();

    message("Done", false);
    mem_peak();
    timer.startTime();

    return 0;
}