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
    log("Pattern length: ", pattern.size());
    log("Pattern: ", (pattern.size() > 100 ? pattern.substr(0, 100) + "..." : pattern));

    std::string pml_filename = args.pattern_filename + ".pml";
    std::string cid_filename = args.pattern_filename + ".cid";
    std::string rev_pml_filename = pml_filename + ".rev";
    std::string rev_cid_filename = cid_filename + ".rev";
    std::ofstream fs_pml(rev_pml_filename);
    std::ofstream fs_cid(rev_cid_filename);
    tbl.query_pml(pattern, fs_pml, fs_cid);

    log("PMLs/CIDs written in reverse, using 'rev' to output final values");
    auto reverse_file = [](const std::string &filename) {
        std::string rev_filename = filename + ".rev";
        std::string cmd = "rev " + rev_filename + " > " + filename + " && rm " + rev_filename;
        log("Command: ", cmd);
        if(system(cmd.c_str()) != 0) {
            error("Failed to reverse file: " + rev_filename);
        }
    };
    reverse_file(pml_filename);
    reverse_file(cid_filename);

    timer.end();
    submessage("Query Complete");
    timer.midTime();
    submessage("PMLs written to: " + pml_filename);
    submessage("CIDs written to: " + cid_filename);

    message("Done", false);
    mem_peak();
    timer.startTime();

    return 0;
}