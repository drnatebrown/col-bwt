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
   \file cid_split.cpp
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

Options::Mode parse_mode(Args &args) {
    Options::Mode mode;
    if (args.split_mode == "default") {
        log("Split mode: Default (find all col runs)");
        mode = Options::Mode::Default;
    } else if (args.split_mode == "tunnels") {
        log("Split mode: Tunnels (find col runs which form BWT tunnels)");
        mode = Options::Mode::Tunneled;
    } else {
        terminal_error("Invalid split mode: " + args.split_mode + ". Must be one of: default, tunnels");
    }
    return mode;
}

Options::Overlap parse_overlap(Args &args) {
    Options::Overlap overlap_mode;
    if (args.overlap == "split") {
        log("Overlap mode: Split (split col runs at overlap)");
        overlap_mode = Options::Overlap::Split;
    } else if (args.overlap == "append") {
        log("Overlap mode: Append (append overlapping col runs after existing)");
        overlap_mode = Options::Overlap::Append;
    } else {
        terminal_error("Invalid overlap mode: " + args.overlap + ". Must be one of: split, append");
    }
    return overlap_mode;
}

int main(int argc, char *const argv[])
{
    Timer timer;
    Args args;
    parseArgs(argc, argv, args);

    if (!args.N) {
        error("Invalid number of COL positions");
        return 1;
    }

    message("Loading BWT table supporting FL mapping: ");
    timer.start();

    FL_table tbl;
    std::string filename_tbl = args.filename + tbl.get_file_extension();
    ifstream fs_tbl(filename_tbl);
    tbl.load(fs_tbl);

    tbl.bwt_stats();

    timer.end();
    submessage("Load Complete");
    timer.startTime();

    timer.mid();
    message("Loading COL Positions: ");

    std::string filename_mums = args.filename + ".mums";
    std::ifstream mum_file(filename_mums, std::ios::binary | std::ios::ate);
    size_t num_mums = mum_file.tellg()/(RW_BYTES*4);

    log("# of COL positions: ", num_mums);
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

    timer.end();
    submessage("Load Complete");
    timer.midTime();

    message("Splitting runs based on COL Positions using FL Table");
    timer.mid();

    Options::Mode mode = parse_mode(args);
    Options::Overlap overlap = parse_overlap(args);
    log("N, # sequences in collection: ", args.N);
    log("Split rate: ", args.split_rate);

    col_split<> split_ds(tbl, mode, overlap);
    split_ds.split(match_lens, match_pos, args.N, args.split_rate);

    timer.end();
    submessage("Splitting Complete");
    timer.midTime();

    message("Serializing COL runs bitvector and IDs");
    timer.mid();

    // split_ds.save(args.filename);

    timer.end();
    submessage("Serialization Complete");
    timer.midTime();

    message("Done", false);
    mem_peak();
    timer.startTime();

    return 0;
}