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
#include <io.hpp>
#include <col_bwt.hpp>
#include <cstddef>
#include <fstream>

#include <sdsl/io.hpp>
#include <sdsl/int_vector.hpp>
#include <malloc_count.h>

void pml_direct_to_file(const std::string &pattern_filename, const std::string &pml_filename, const std::string &cid_filename, col_pml &tbl)
{
    std::string rev_pml_filename = pml_filename + ".rev";
    std::string rev_cid_filename = cid_filename + ".rev";
    std::ofstream fs_pml_rev(rev_pml_filename);
    std::ofstream fs_cid_rev(rev_cid_filename);

    PatternProcessor patterns(pattern_filename);
    while (patterns.read()) {
        const auto& [pattern_s, m] = patterns.get_seq();
        std::string header_rev = '>' + patterns.get_id() + " \n";
        std::reverse(header_rev.begin(), header_rev.end());

        fs_pml_rev << header_rev;
        fs_cid_rev << header_rev;
        tbl.query_pml(pattern_s, m, fs_pml_rev, fs_cid_rev);
        fs_pml_rev << "\n";
        fs_cid_rev << "\n";
    }

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
}

void pml_to_vec(const std::string &pattern_filename, const std::string &pml_filename, const std::string &cid_filename, col_pml &tbl)
{
    std::ofstream fs_pml(pml_filename);
    std::ostream_iterator<size_t> pml_iter(fs_pml, " ");
    
    std::ofstream fs_cid(cid_filename);
    std::ostream_iterator<size_t> cid_iter(fs_cid, " ");

    PatternProcessor patterns(pattern_filename);
    while (patterns.read()) {
        const auto& [pattern_s, m] = patterns.get_seq();
        const auto& [pmls, cids] = tbl.query_pml(pattern_s, m);

        // TODO - why does spumoni have this extra space?
        fs_pml << '>' << patterns.get_id() << " \n";
        std::copy(pmls.begin(), pmls.end(), pml_iter);
        fs_pml << "\n";

        fs_cid << '>' << patterns.get_id() << " \n";
        std::copy(cids.begin(), cids.end(), cid_iter);
        fs_cid << "\n";
    }
    
    fs_pml.close();
    fs_cid.close();
}

int main(int argc, char *const argv[])
{
    Timer timer;
    Args args;
    parseArgs(argc, argv, args);

    if (args.filename.empty())
    {
        error("BWT table not provided");
    }
    if (args.pattern_filename.empty())
    {
        error("Pattern file not provided");
    }

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
    message("Computing PML Queries: ");
    
    std::string pml_filename = args.pattern_filename + ".pml";
    std::string cid_filename = args.pattern_filename + ".cid";
    if (args.long_pattern) {
        pml_direct_to_file(args.pattern_filename, pml_filename, cid_filename, tbl);
    } else {
        pml_to_vec(args.pattern_filename, pml_filename, cid_filename, tbl);
    }

    // message("Reading pattern file: ");

    // std::ifstream fs_pattern(args.pattern_filename);
    // std::string pattern;
    // std::string line;
    // std::ostringstream sequence_stream;
    // while (std::getline(fs_pattern, line))
    // {
    //     if (line.empty() || line[0] == '>')
    //         continue;
    //     sequence_stream << line;
    // }
    // pattern = sequence_stream.str();

    // timer.end();
    // submessage("Read Complete");
    // timer.midTime();

    // timer.mid();
    // message("Computing PMLs Query: ");
    // log("Pattern length: ", pattern.size());
    // log("Pattern: ", (pattern.size() > 100 ? pattern.substr(0, 100) + "..." : pattern));

    // std::string pml_filename = args.pattern_filename + ".pml";
    // std::string cid_filename = args.pattern_filename + ".cid";
    // std::string rev_pml_filename = pml_filename + ".rev";
    // std::string rev_cid_filename = cid_filename + ".rev";
    // std::ofstream fs_pml(rev_pml_filename);
    // std::ofstream fs_cid(rev_cid_filename);
    // tbl.query_pml(pattern, fs_pml, fs_cid);

    // log("PMLs/CIDs written in reverse, using 'rev' to output final values");
    // auto reverse_file = [](const std::string &filename) {
    //     std::string rev_filename = filename + ".rev";
    //     std::string cmd = "rev " + rev_filename + " > " + filename + " && rm " + rev_filename;
    //     log("Command: ", cmd);
    //     if(system(cmd.c_str()) != 0) {
    //         error("Failed to reverse file: " + rev_filename);
    //     }
    // };
    // reverse_file(pml_filename);
    // reverse_file(cid_filename);

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