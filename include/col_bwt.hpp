/* col_bwt - Uncompressed version of OptBWTR (LF table) with co-lin ids
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
   \file col_bwt.hpp
   \brief Uncompressed version of OptBWTR (LF table) with co-lin ids
   \author Nathaniel Brown
   \date 07/08/2024
*/

#ifndef COL_BWT_HPP
#define COL_BWT_HPP

#include "LF_table.hpp"
#include "common/common.hpp"

#include <common.hpp>
#include <iostream>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>
#include <sdsl/sd_vector.hpp>
#include <sys/types.h>

using namespace std;

// Row of the LF table
class col_row : public LF_row
{
public:
    ulint col_id : ID_BITS;

    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name ="")
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        size_t temp = col_id;
        out.write((char *)&temp, ID_BYTES);
        written_bytes += ID_BYTES;

        LF_row::serialize(out, child, "LF_row");

        return written_bytes;
    }

    void load(std::istream &in)
    {
        size_t temp;
        in.read((char *)&temp, BWT_BYTES);
        col_id = temp;

        LF_row::load(in);
    }
};

class col_bwt : public LF_table<col_row>
{
public:
    col_bwt() {}

    col_bwt(std::ifstream &heads, std::ifstream &lengths, std::ifstream &col_ids, sdsl::sd_vector<> splits)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);
        col_ids.clear();
        col_ids.seekg(0);
        
        LF_runs = vector<col_row>();
        vector<vector<size_t>> L_block_indices = vector<vector<size_t>>(ALPHABET_SIZE);
        
        size_t s_set_bits = 0;
        {
            sd_rank s_rank(&splits);
            s_set_bits = s_rank(splits.size());
            assert(s_set_bits > 1);

            col_ids.seekg(0, std::ios::end);
            size_t col_ids_size = col_ids.tellg();
            col_ids_size /= ID_BYTES;
            col_ids.seekg(0, std::ios::beg);
            assert(col_ids_size == s_set_bits);
        }
        sd_select s_select(&splits);

        status("Constructing Col BWT table");
        char c;
        size_t i = 0;
        r = 0;
        bwt_r = 0;
        n = 0;
        size_t s = 0;
        size_t s_curr = s_select(s + 1);
        while ((c = heads.get()) != EOF)
        {
            size_t length = 0;
            lengths.read((char *)&length, 5);
            if (c <= TERMINATOR) c = TERMINATOR;

            while (s < s_set_bits && s_curr < n + length)
            {
                size_t curr_id = 0;
                L_block_indices[c].push_back(i++);

                // If the split corresponds to an existing run-head, write the run with no-id first
                if (s_curr > n) {
                    LF_runs.push_back({c, n, 0, 0, curr_id});
                    col_ids.read((char *)&curr_id, ID_BYTES);

                    ulint delta = s_curr - n;
                    n += delta;
                    length -= delta;
                    ++s;
                    s_curr = (s < s_set_bits) ? s_select(s + 1) : 0;
                }
                // otherwise the split corresponds to this run, fetch its id first
                else {
                    col_ids.read((char *)&curr_id, ID_BYTES);
                    LF_runs.push_back({c, n, 0, 0, curr_id});

                    ++s;
                    s_curr = (s < s_set_bits) ? s_select(s + 1) : 0;
                    // The length of this run is determined by the next split
                    if (s < s_set_bits && s_curr < n + length) {
                        ulint delta = s_curr - n;
                        n += delta;
                        length -= delta;
                    }
                    // Else its the rest of this run
                    else {
                        n += length;
                        length = 0;
                    }
                }
            }

            // Regular BWT run
            if (length > 0) {
                LF_runs.push_back({c, n, 0, 0, 0});
                L_block_indices[c].push_back(i++);
                n+=length;
            }

            ++bwt_r;
        }
        r = LF_runs.size();
        status();

        status("Computing values of Col BWT table");
        compute_table(L_block_indices);
        status();

        #ifdef PRINT_STATS
        stat("Text runs", runs());
        stat("Text length", size());
        #endif
    }

    col_bwt(std::ifstream &bwt, std::ifstream &col_ids, sdsl::sd_vector<> splits)
    {
        bwt.clear();
        bwt.seekg(0);
        col_ids.clear();
        col_ids.seekg(0);
        
        LF_runs = vector<col_row>();
        vector<vector<size_t>> L_block_indices = vector<vector<size_t>>(ALPHABET_SIZE);

        size_t s_set_bits = 0;
        {
            sd_rank s_rank(&splits);
            s_set_bits = s_rank(splits.size());
            assert(s_set_bits > 1);

            col_ids.seekg(0, std::ios::end);
            size_t col_ids_size = col_ids.tellg()/ID_BYTES;
            col_ids.seekg(0, std::ios::beg);
            assert(col_ids_size == s_set_bits);
        }
        sd_select s_select(&splits);
        
        status("Constructing Col BWT table");
        char last_c;
        char c;
        size_t length = 0;
        size_t i = 0;
        r = 0;
        bwt_r = 1;
        n = 0;
        size_t s = 0;
        size_t s_curr = s_select(s + 1);
        size_t curr_id = 0;
        while ((c = bwt.get()) != EOF)
        {
            if (c <= TERMINATOR) c = TERMINATOR;
            
            if (i != 0 && (c != last_c || i == s_curr))
            {
                LF_runs.push_back({last_c, n, 0, 0, curr_id});
                L_block_indices[last_c].push_back(i++);

                n+=length;
                length = 0;
                
                if (i == s_curr) {
                    col_ids.read((char *)&curr_id, ID_BYTES);
                    ++s;
                    if (s < s_set_bits) s_curr = s_select(s + 1);
                }
                else {
                    curr_id = 0;
                }

                if (c != last_c) {
                    ++bwt_r;
                }
            }
            ++length;
            last_c = c;
        }
        // Step for final character
        LF_runs.push_back({last_c, length, 0, 0, curr_id});
        L_block_indices[last_c].push_back(i++);  
        n+=length;

        r = LF_runs.size();
        status();

        status("Computing values of Col BWT table");
        compute_table(L_block_indices);
        status();

        #ifdef PRINT_STATS
        stat("Text runs", runs());
        stat("Text length", size());
        #endif
    }

    ulint bwt_runs() const
    {
        return bwt_r;
    }

    std::string get_file_extension() const
    {
        return ".col_bwt";
    }

    void bwt_stats()
    {
        log("Number of Col equal-letter runs: r = ", r);
        log("Number of BWT equal-letter runs: r = ", bwt_r);
        log("Length of complete BWT: n = ", n);
        log("Rate n/r = ", double(n) / r);
        log("log2(r) = ", log2(double(r)));
        log("log2(n/r) = ", log2(double(n) / r));
    }

    void mem_stats()
    {
        log("Memory (bytes):");
        log("   Col BWT: ", r + 4*r*BWT_BYTES + ID_BYTES);
        log("           Chars: ", r);
        log("            Ints: ", r*BWT_BYTES);
        log("             IDs: ", r*ID_BYTES);
    }

    /* serialize to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name ="")
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        out.write((char *)&bwt_r, sizeof(bwt_r));
        written_bytes += sizeof(bwt_r);

        written_bytes += LF_table::serialize(out);

        return written_bytes;
    }

    /* load from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        size_t size;

        in.read((char *)&bwt_r, sizeof(bwt_r));

        LF_table::load(in);
    }

protected:
    ulint bwt_r;
};

#endif /* end of include guard: _COL_BWT_HH */