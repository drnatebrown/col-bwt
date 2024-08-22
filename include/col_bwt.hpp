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
#include <cstddef>
#include <thread>
#include <iostream>
#include <functional>
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

    col_row() {}

    col_row(uchar c, ulint ix, ulint in, ulint o, ulint id)
        : LF_row(c, ix, in, o), col_id(id) {}

    size_t serialize(std::ostream &out)
    {
        size_t written_bytes = 0;

        size_t temp = col_id;
        out.write((char *)&temp, ID_BYTES);
        written_bytes += ID_BYTES;

        LF_row::serialize(out);

        return written_bytes;
    }

    void load(std::istream &in)
    {
        size_t temp;
        in.read((char *)&temp, BWT_BYTES);
        col_id = temp;

        LF_row::load(in);
    }
} __attribute__((packed));

// Row of the LF table using thresholds
class col_thr : public col_row
{
public:
    ulint threshold : BWT_BITS;

    col_thr() {}

    col_thr(uchar c, ulint ix, ulint in, ulint o, ulint id, ulint t)
        : col_row(c, ix, in, o, id), threshold(t) {}

    col_thr(uchar c, ulint ix, ulint in, ulint o, ulint id)
        : col_row(c, ix, in, o, id), threshold() {}

    size_t serialize(std::ostream &out)
    {
        size_t written_bytes = 0;

        size_t temp = threshold;
        out.write((char *)&temp, BWT_BYTES);
        written_bytes += BWT_BYTES;

        col_row::serialize(out);

        return written_bytes;
    }

    void load(std::istream &in)
    {
        size_t temp;
        in.read((char *)&temp, BWT_BYTES);
        threshold = temp;

        col_row::load(in);
    }
} __attribute__((packed));

// Must use col_row or an inherited class of col_row
template <class row_t = col_row>
class col_bwt : public LF_table<row_t>
{
public:
    col_bwt() {}

    col_bwt(std::ifstream &heads, std::ifstream &lengths, std::ifstream &col_ids, sdsl::sd_vector<> &splits)
    {
        #ifdef PRINT_STATS
        ulint col_chars = 0;
        ulint col_runs = 0;
        #endif

        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);
        col_ids.clear();
        col_ids.seekg(0);
        
        this->LF_runs = vector<row_t>();
        vector<vector<size_t>> L_block_indices = vector<vector<size_t>>(ALPHABET_SIZE);
        
        size_t s_set_bits = 0;
        {
            sd_rank s_rank(&splits);
            s_set_bits = s_rank(splits.size());
            assert(s_set_bits > 1);

            #ifdef DEBUG
            col_ids.seekg(0, std::ios::end);
            size_t col_ids_size = col_ids.tellg();
            col_ids_size /= ID_BYTES;
            col_ids.clear();
            col_ids.seekg(0, std::ios::beg);
            assert(col_ids_size == s_set_bits);
            #endif
        }
        sd_select s_select(&splits);

        status("Constructing Col BWT table");
        char c;
        size_t i = 0;
        this->r = 0;
        bwt_r = 0;
        this->n = 0;
        size_t s = 0;
        size_t s_curr = s_select(s + 1);
        size_t curr_id = 0;
        while ((c = heads.get()) != EOF)
        {
            size_t length = 0;
            lengths.read((char *)&length, RW_BYTES);
            if (c <= TERMINATOR) c = TERMINATOR;
            #ifdef DNA_ALPHABET
            c = charToBits[c];
            #endif

            size_t run_end = this->n + length;
            if (s_curr == this->n) {
                col_ids.read((char *)&curr_id, ID_BYTES);
                ++s;
                s_curr = (s < s_set_bits) ? s_select(s + 1) : 0;
            }

            while (s < s_set_bits && s_curr < run_end) {
                L_block_indices[c].push_back(i++);
                this->LF_runs.push_back(row_t(c, this->n, 0, 0, curr_id));
                size_t delta = s_curr - this->n;
                this->n += delta;
                length -= delta;

                #ifdef PRINT_STATS
                if (curr_id > 0) {
                    col_chars += delta;
                    ++col_runs;
                }
                #endif

                ++s;
                s_curr = (s < s_set_bits) ? s_select(s + 1) : 0;
                col_ids.read((char *)&curr_id, ID_BYTES);
            }

            if (length > 0) {
                L_block_indices[c].push_back(i++);
                this->LF_runs.push_back(row_t(c, this->n, 0, 0, curr_id));
                this->n += length;

                #ifdef PRINT_STATS
                if (curr_id > 0) {
                    col_chars += length;
                    ++col_runs;
                }
                #endif
            }
            ++bwt_r;
        }
        this->r = this->LF_runs.size();
        status();

        status("Computing values of Col BWT table");
        this->compute_table(L_block_indices);
        status();

        #ifdef PRINT_STATS
        stat("Col-BWT runs", this->runs());
        stat("Col runs", col_runs);
        stat("Col chars", col_chars);
        stat("BWT runs", bwt_runs());
        stat("Text length", this->size());
        #endif
    }

    col_bwt(std::ifstream &bwt, std::ifstream &col_ids, sdsl::sd_vector<> splits)
    {
        #ifdef PRINT_STATS
        ulint col_chars = 0;
        ulint col_runs = 0;
        #endif

        bwt.clear();
        bwt.seekg(0);
        col_ids.clear();
        col_ids.seekg(0);
        
        this->LF_runs = vector<row_t>();
        vector<vector<size_t>> L_block_indices = vector<vector<size_t>>(ALPHABET_SIZE);

        size_t s_set_bits = 0;
        {
            sd_rank s_rank(&splits);
            s_set_bits = s_rank(splits.size());
            assert(s_set_bits > 1);

            #ifdef DEBUG
            col_ids.seekg(0, std::ios::end);
            size_t col_ids_size = col_ids.tellg()/ID_BYTES;
            col_ids.seekg(0, std::ios::beg);
            assert(col_ids_size == s_set_bits);
            #endif
        }
        sd_select s_select(&splits);

        status("Constructing Col BWT table");
        char last_c;
        char c;
        size_t length = 0;
        size_t i = 0;
        this->r = 0;
        bwt_r = 1;
        this->n = 0;
        size_t s = 0;
        size_t s_curr = s_select(s + 1);
        size_t curr_id = 0;
        while ((c = bwt.get()) != EOF)
        {
            if (c <= TERMINATOR) c = TERMINATOR;
            #ifdef DNA_ALPHABET
            c = charToBits[c];
            #endif
            
            if (i != 0 && (c != last_c || i == s_curr))
            {
                this->LF_runs.push_back(row_t(last_c, this->n, 0, 0, curr_id));
                L_block_indices[last_c].push_back(i++);

                #ifdef PRINT_STATS
                if (curr_id > 0) {
                    col_chars += length;
                    ++col_runs;
                }
                #endif

                this->n+=length;
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
        this->LF_runs.push_back(row_t(last_c, length, 0, 0, curr_id));
        L_block_indices[last_c].push_back(i++);  
        this->n += length;

        this->r = this->LF_runs.size();
        status();

        status("Computing values of Col BWT table");
        this->compute_table(L_block_indices);
        status();

        #ifdef PRINT_STATS
        stat("Col-BWT runs", this->runs());
        stat("Col runs", col_runs);
        stat("BWT runs", bwt_runs());
        stat("Text length", this->size());
        #endif
    }

    ulint bwt_runs() const
    {
        return bwt_r;
    }

    void bwt_stats()
    {
        log("Number of Col equal-letter runs: r = ", this->r);
        log("Number of BWT equal-letter runs: bwt_r = ", bwt_r);
        log("Length of complete BWT: n = ", this->n);
        log("Rate n/r = ", double(this->n) / this->r);
        log("log2(r) = ", log2(double(this->r)));
        log("log2(n/r) = ", log2(double(this->n) / this->r));
    }

    void mem_stats()
    {
        log("Memory (bytes):");
        log("   Col BWT: ", bits_to_bytes(this->r*ALPHABET_BITS + 3*this->r*BWT_BITS + this->r*ID_BITS));
        log("           Chars: ", bits_to_bytes(this->r*ALPHABET_BITS));
        log("            Ints: ", bits_to_bytes(this->r*BWT_BITS));
        log("             IDs: ", bits_to_bytes(this->r*ID_BITS));
    }

    std::string get_file_extension() const
    {
        return ".col_bwt";
    }

    /* serialize to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name ="")
    {
        size_t written_bytes = 0;

        out.write((char *)&bwt_r, sizeof(bwt_r));
        written_bytes += sizeof(bwt_r);

        written_bytes += LF_table<row_t>::serialize(out);

        return written_bytes;
    }

    /* load from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        in.read((char *)&bwt_r, sizeof(bwt_r));

        LF_table<row_t>::load(in);
    }

protected:
    ulint bwt_r;
};

class col_pml : public col_bwt<col_thr>
{
public:
    col_pml() {}

    col_pml(std::ifstream &heads, std::ifstream &lengths, std::ifstream &col_ids, std::ifstream &thresholds, sdsl::sd_vector<> splits)
        : col_bwt<col_thr>(heads, lengths, col_ids, splits)
    {
        read_thresholds(thresholds);
    }

    col_pml(std::ifstream &bwt, std::ifstream &col_ids, std::ifstream &thresholds, sdsl::sd_vector<> splits)
        : col_bwt<col_thr>(bwt, col_ids, splits)
    {
        read_thresholds(thresholds);
    }

    std::pair<std::vector<ulint>, std::vector<ulint>> query_pml(const std::string &pattern)
    {
        size_t m = pattern.size();
        return _query_pml(pattern.data(), m);
    }

    std::pair<std::vector<ulint>, std::vector<ulint>> query_pml(const char *pattern, const size_t m)
    {
        return _query_pml(pattern, m);
    }

    void query_pml(const std::string &pattern, std::ofstream &lens_rev, std::ofstream &col_ids_rev)
    {
        size_t m = pattern.size();
        _query_pml(pattern.data(), m, lens_rev, col_ids_rev);
    }

    void query_pml(const char *pattern, const size_t m, std::ofstream &lens_rev, std::ofstream &col_ids_rev)
    {
        _query_pml(pattern, m, lens_rev, col_ids_rev);
    }

    void mem_stats()
    {
        log("Memory (bytes):");
        log("   Col BWT: ", bits_to_bytes(this->r*ALPHABET_BITS + 4*this->r*BWT_BITS + this->r*ID_BITS));
        log("           Chars: ", bits_to_bytes(this->r*ALPHABET_BITS));
        log("            Ints: ", bits_to_bytes(this->r*BWT_BITS));
        log("             IDs: ", bits_to_bytes(this->r*ID_BITS));
    }

    std::string get_file_extension() const
    {
        return ".col_pml";
    }

protected:
    void read_thresholds(std::ifstream &thresholds)
    {
        status("Reading thresholds");
        thresholds.clear();
        thresholds.seekg(0);

        size_t threshold = 0;
        size_t i = 0;
        while (thresholds.read((char *)&threshold, RW_BYTES))
        {
            uchar last = get_char(i);
            do {
                this->LF_runs[i++].threshold = threshold;
            } while (i < this->r && get_char(i) == last);
        }
        assert(i == this->bwt_r);
        status();
    }

    // Returns vectors of PML lengths and column ids
    std::pair<vector<ulint>, vector<ulint>> _query_pml(const char *pattern, const size_t m)
    {
        std::vector<ulint> pml_lengths(m);
        std::vector<ulint> col_stats(m);
        auto store_to_vector = [&](const ulint length, const ulint col_id, const ulint idx) {
            pml_lengths[idx] = length;
            col_stats[idx] = col_id;
        };

        _query_pml(pattern, m, store_to_vector);

        return std::make_pair(pml_lengths, col_stats);
    }

    // Writes to file in reverse
    void _query_pml(const char *pattern, const size_t m, std::ofstream &lens_rev, std::ofstream &col_ids_rev) 
    {
        lens_rev.clear();
        lens_rev.seekp(0);
        col_ids_rev.clear();
        col_ids_rev.seekp(0);

        auto store_to_file_rev = [&](const ulint length, const ulint col_id, const ulint idx) {
            std::string len_str = std::to_string(length);
            std::string col_id_str = std::to_string(col_id);
            std::reverse(len_str.begin(), len_str.end());
            std::reverse(col_id_str.begin(), col_id_str.end());

            lens_rev << " " << len_str;
            col_ids_rev << " " << col_id_str;
        };

        _query_pml(pattern, m, store_to_file_rev);

        lens_rev.close();
        col_ids_rev.close();
    }

    void _query_pml(const char *pattern, const size_t m, std::function<void(const ulint, const ulint, const ulint)> store_stats)
    {
        std::vector<ulint> pml_lengths(m);
        std::vector<ulint> col_stats(m);

        ulint pos = this->n - 1;
        ulint interval = r - 1;
        ulint offset = get_length(interval) - 1;

        ulint length = 0;
        ulint col_id = 0;

        for (size_t i = 0; i < m; ++i)
        {
            uchar c = pattern[m - i - 1];
            col_id = LF_runs[interval].col_id;

            // Case 1: c is the same as the current character, extend match
            if (get_char(interval) == c) {
                ++length;
            }
            // Case 2: mismatch, reset PML to 0 and use thresholds to reorient in BWT
            else {
                length = 0;
                threshold_step(interval, offset, pos, c);
            }

            store_stats(length, col_id, m - i - 1);

            std::tie(interval, offset, pos) = LF_idx(interval, offset); // Next step
        }
    }

    void threshold_step(ulint &interval, ulint &offset, const ulint pos, const uchar c)
    {
        ulint new_interval = interval;
        ulint new_offset = offset;
        ulint thr = this->n;

        #ifdef MULTI_THREAD
        std::optional<std::pair<ulint, ulint>> pred_result, succ_result;
        std::thread pred_thread([&]() {
            pred_result = pred_char(interval, c);
        });
        std::thread succ_thread([&]() {
            succ_result = succ_char(interval, c);
        });
        pred_thread.join();
        succ_thread.join();
        #else
        std::optional<std::pair<ulint, ulint>> succ_result = succ_char(interval, c);
        #endif

        // Has successor, get threshold
        if (succ_result.has_value()) {
            auto [succ_interval, succ_offset] = succ_result.value();
            thr = LF_runs[succ_interval].threshold;
            new_interval = succ_interval;
            new_offset = succ_offset;
        }

        // Predecessor is better match, or no successor
        if (pos < thr) {
            #ifndef MULTI_THREAD
            std::optional<std::pair<ulint, ulint>> pred_result = pred_char(interval, c);
            #endif

            if (pred_result.has_value()) {
                auto[pred_interval, pred_offset] = pred_result.value();
                new_interval = pred_interval;
                new_offset = pred_offset;
            }
        }
        
        interval = new_interval;
        offset = new_offset;
    }
};
#endif /* end of include guard: _COL_BWT_HH */