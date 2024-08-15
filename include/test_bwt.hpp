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

#ifndef TEST_BWT_HPP
#define TEST_BWT_HPP

#include "common/common.hpp"

#include <common.hpp>
#include <fstream>
#include <thread>
#include <iostream>
#include <functional>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>
#include <sdsl/sd_vector.hpp>
#include <sys/types.h>

using namespace std;

// Row of the LF table using thresholds
struct __attribute__((packed)) test_row {
    uchar character : ALPHABET_BITS;
    ulint idx : BWT_BITS;
    ulint interval : BWT_BITS;
    ulint offset : BWT_BITS;
    ulint col_id : ID_BITS;
    ulint threshold : BWT_BITS;

    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name ="")
    {
        size_t total_size = bits_to_bytes(ALPHABET_BITS + 4*BWT_BITS + ID_BITS);

        // Allocate a byte array
        char *buffer = new char[total_size];
        
        // Copy each member into the byte array
        size_t cursor = 0;
        size_t temp = character;
        memcpy(buffer + cursor, &temp, sizeof(uchar));
        cursor += sizeof(uchar);

        temp = idx;
        memcpy(buffer + cursor, &temp, BWT_BYTES);
        cursor += BWT_BYTES;

        temp = interval;
        memcpy(buffer + cursor, &temp, BWT_BYTES);
        cursor += BWT_BYTES;

        temp = offset;
        memcpy(buffer + cursor, &temp, BWT_BYTES);
        cursor += BWT_BYTES;

        temp = col_id;
        memcpy(buffer + cursor, &temp, ID_BYTES);
        cursor += ID_BYTES;

        temp = threshold;
        memcpy(buffer + cursor, &temp, BWT_BYTES);

        // Write the byte array to the output stream
        out.write(buffer, total_size);
        
        // Clean up
        delete[] buffer;
        
        return total_size;
    }
    
    // {
    //     sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    //     size_t written_bytes = 0;

    //     size_t temp = character;
    //     out.write((char *)&temp, sizeof(char));
    //     written_bytes += sizeof(char);

    //     temp = idx;
    //     out.write((char *)&temp, BWT_BYTES);
    //     written_bytes += BWT_BYTES;

    //     temp = interval;
    //     out.write((char *)&temp, BWT_BYTES);
    //     written_bytes += BWT_BYTES;

    //     temp = offset;
    //     out.write((char *)&temp, BWT_BYTES);
    //     written_bytes += BWT_BYTES;

    //     temp = col_id;
    //     out.write((char *)&temp, ID_BYTES);
    //     written_bytes += ID_BYTES;

    //     temp = threshold;
    //     out.write((char *)&temp, BWT_BYTES);
    //     written_bytes += BWT_BYTES;

    //     return written_bytes;
    // }

    void load(std::istream &in)
    {
        size_t total_size = sizeof(uchar) + 4 * BWT_BYTES + ID_BYTES;
        char buffer[total_size];

        in.read(buffer, total_size);

        size_t cursor = 0;
        size_t temp;

        memcpy(&temp, buffer + cursor, sizeof(uchar));
        character = temp;
        cursor += sizeof(uchar);

        memcpy(&temp, buffer + cursor, BWT_BYTES);
        idx = temp;
        cursor += BWT_BYTES;

        memcpy(&temp, buffer + cursor, BWT_BYTES);
        interval = temp;
        cursor += BWT_BYTES;

        memcpy(&temp, buffer + cursor, BWT_BYTES);
        offset = temp;
        cursor += BWT_BYTES;

        memcpy(&temp, buffer + cursor, ID_BYTES);
        col_id = temp;
        cursor += ID_BYTES;

        memcpy(&temp, buffer + cursor, BWT_BYTES);
        threshold = temp;

        // size_t temp;
        // in.read((char *)&temp, sizeof(char));
        // character = temp;

        // in.read((char *)&temp, BWT_BYTES);
        // idx = temp;

        // in.read((char *)&temp, BWT_BYTES);
        // interval = temp;

        // in.read((char *)&temp, BWT_BYTES);
        // offset = temp;

        // in.read((char *)&temp, ID_BYTES);
        // col_id = temp;

        // in.read((char *)&temp, BWT_BYTES);
        // threshold = temp;
    }
};

template <class row_t = test_row>
class test_bwt
{
public:
    test_bwt() {}

    test_bwt(std::ifstream &heads, std::ifstream &lengths, std::ifstream &col_ids, std::ifstream &thresholds,sdsl::sd_vector<> splits)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);
        col_ids.clear();
        col_ids.seekg(0);
        thresholds.clear();
        thresholds.seekg(0);
        
        LF_runs = vector<row_t>();
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
            lengths.read((char *)&length, RW_BYTES);
            size_t threshold = 0;
            thresholds.read((char *)&threshold, RW_BYTES);
            if (c <= TERMINATOR) c = TERMINATOR;
            #ifdef DNA_ALPHABET
            c = charToBits[c];
            #endif

            while (s < s_set_bits && s_curr < n + length)
            {
                size_t curr_id = 0;
                L_block_indices[c].push_back(i++);

                // If the split corresponds to an existing run-head, write the run with no-id first
                if (s_curr >n) {
                    LF_runs.push_back({c, n, 0, 0, curr_id, threshold});
                    col_ids.read((char *)&curr_id, ID_BYTES);

                    ulint delta = s_curr -n;
                    n += delta;
                    length -= delta;
                    ++s;
                    s_curr = (s < s_set_bits) ? s_select(s + 1) : 0;
                }
                // otherwise the split corresponds to this run, fetch its id first
                else {
                    col_ids.read((char *)&curr_id, ID_BYTES);
                    LF_runs.push_back({c, n, 0, 0, curr_id, threshold});

                    ++s;
                    s_curr = (s < s_set_bits) ? s_select(s + 1) : 0;
                    // The length of this run is determined by the next split
                    if (s < s_set_bits && s_curr <n + length) {
                        ulint delta = s_curr -n;
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
                LF_runs.push_back({c, n, 0, 0, 0, threshold});
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
        stat("Col runs", runs());
        stat("BWT runs", bwt_runs());
        stat("Text length", size());
        #endif
    }

    test_bwt(std::ifstream &bwt, std::ifstream &col_ids, std::ifstream &thresholds, sdsl::sd_vector<> splits)
    {
        bwt.clear();
        bwt.seekg(0);
        col_ids.clear();
        col_ids.seekg(0);
        thresholds.clear();
        thresholds.seekg(0);
        
        LF_runs = vector<row_t>();
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
        size_t threshold = 0;
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
            #ifdef DNA_ALPHABET
            c = charToBits[c];
            #endif

            if (i != 0 && (c != last_c || i == s_curr))
            {
                if (c != last_c) {
                    ++bwt_r;
                    thresholds.read((char *)&threshold, RW_BYTES);
                }

                LF_runs.push_back({last_c, n, 0, 0, curr_id, threshold});
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
            }
            ++length;
            last_c = c;
        }
        thresholds.read((char *)&threshold, RW_BYTES);
        // Step for final character
        LF_runs.push_back({last_c, n, 0, 0, curr_id, threshold});
        L_block_indices[last_c].push_back(i++);  
        n += length;

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
    const row_t get(size_t i)
    {
        assert(i < LF_runs.size());
        return LF_runs[i];
    }

    uchar get_char(ulint i)
    {
        #ifdef DNA_ALPHABET
        return bitsToChar[get(i).character];
        #else
        return get(i).character;
        #endif
    }

    uchar get_char_bits(ulint i)
    {
        return get(i).character;
    }

    // TODO specilization for template with lengths, not idx
    ulint get_length(ulint i)
    {
        return (i == r - 1) ? (n - get_idx(i)) : (get_idx(i + 1) - get_idx(i));
    }

    ulint get_idx(ulint i)
    {
        return get(i).idx;
    }

    ulint to_idx(ulint interval, ulint offset)
    {
        return get_idx(interval) + offset;
    }

    ulint size()
    {
        return n;
    }

    ulint runs()
    {
        return r;
    }

    void invert(std::string outfile) 
    {
        std::ofstream out(outfile);

        ulint interval = 0;
        ulint offset = 0;

        uchar c;
        while((c = get_char(interval)) > TERMINATOR) 
        {
            out << c;
            std::pair<ulint, ulint> pos = LF(interval, offset);
            interval = pos.first;
            offset = pos.second;
        }
    }

    /*
     * \param Run position (RLE intervals)
     * \param Current character offset in block
     * \return block position and offset of preceding character
     */
    std::pair<ulint, ulint> LF(ulint run, ulint offset)
    {
        ulint next_interval = LF_runs[run].interval;
        ulint next_offset = LF_runs[run].offset + offset;

        while (next_offset >= get_length(next_interval)) 
        {
            next_offset -= get_length(next_interval++);
        }

        return std::make_pair(next_interval, next_offset);
    }

    std::tuple<ulint, ulint, ulint> LF_idx(ulint run, ulint offset)
    {
        auto [next_interval, next_offset] = LF(run, offset);
        return std::make_tuple(next_interval, next_offset, to_idx(next_interval, next_offset));
    }

    /* Returns row of largest idx before or at position run with character c */
    std::optional<std::pair<ulint, ulint>> pred_char(ulint run, uchar c)
    {
        #ifdef DNA_ALPHABET
        c = charToBits[c];
        #endif
        while (get_char_bits(run) != c) 
        {
            if (run == 0) return std::nullopt;
            --run;
        }

        return std::make_pair(run, get_length(run) - 1);
    }

    /* Returns row of smallest idx after or at position run with character c */
    std::optional<std::pair<ulint, ulint>> succ_char(ulint run, uchar c)
    {
        #ifdef DNA_ALPHABET
        c = charToBits[c];
        #endif
        while (get_char_bits(run) != c)  
        {
            if (run == r - 1) return std::nullopt;
            ++run;
        }

        return std::make_pair(run, 0);
    }

    void bwt_stats()
    {
        log("Number of Col equal-letter runs: r = ", r);
        log("Number of BWT equal-letter runs: bwt_r = ", bwt_r);
        log("Length of complete BWT: n = ", n);
        log("Rate n/r = ", double(n) / r);
        log("log2(r) = ", log2(double(r)));
        log("log2(n/r) = ", log2(double(n) / r));
    }

    void mem_stats()
    {
        log("Memory (bytes):");
        log("   Col BWT: ", bits_to_bytes(r*ALPHABET_BITS + 4*r*BWT_BITS + r*ID_BITS));
        log("           Chars: ", bits_to_bytes(r*ALPHABET_BITS));
        log("            Ints: ", bits_to_bytes(r*BWT_BITS));
        log("             IDs: ", bits_to_bytes(r*ID_BITS));
    }

    std::string get_file_extension() const
    {
        return ".test_bwt";
    }

    /* serialize to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name ="")
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        out.write((char *)&n, sizeof(n));
        written_bytes += sizeof(n);

        out.write((char *)&r, sizeof(r));
        written_bytes += sizeof(r);

        out.write((char *)&bwt_r, sizeof(bwt_r));
        written_bytes += sizeof(bwt_r);

        size_t size = LF_runs.size();
        out.write((char *)&size, sizeof(size));
        written_bytes += sizeof(size);

        const char* LF_runs_mem = reinterpret_cast<const char*>(LF_runs.data());
        size_t LF_runs_bytes = LF_runs.size() * sizeof(row_t);
        out.write(LF_runs_mem, LF_runs_bytes);
        written_bytes += LF_runs_bytes;

        return written_bytes;
    }

    /* load from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        size_t size;

        in.read((char *)&n, sizeof(n));
        in.read((char *)&r, sizeof(r));
        in.read((char *)&bwt_r, sizeof(bwt_r));

        in.read((char *)&size, sizeof(size));
        LF_runs = std::vector<row_t>(size);

        char* LF_runs_mem = reinterpret_cast<char*>(LF_runs.data());
        size_t LF_runs_bytes = LF_runs.size() * sizeof(row_t);
        in.read(LF_runs_mem, LF_runs_bytes);
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

protected:
    ulint n; // Length of BWT
    ulint r; // Runs of BWT
    ulint bwt_r;

    vector<row_t> LF_runs;

    void compute_table(vector<vector<ulint>> L_block_indices) {\
        ulint curr_L_num = 0;
        ulint L_seen = 0;
        ulint F_seen = 0;
        for(size_t i = 0; i < L_block_indices.size(); ++i) 
        {
            for(size_t j = 0; j < L_block_indices[i].size(); ++j) 
            {
                ulint pos = L_block_indices[i][j];

                LF_runs[pos].interval = curr_L_num;
                LF_runs[pos].offset = F_seen - L_seen;

                F_seen += get_length(pos);

                while (curr_L_num < r && F_seen >= L_seen + get_length(curr_L_num)) 
                {
                    L_seen += get_length(curr_L_num);
                    ++curr_L_num;
                }
            }
        }
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
            lens_rev << ((idx != m - 1) ? "," : "") << len_str;
            col_ids_rev << ((idx != m - 1) ? "," : "") << col_id_str;
        };

        _query_pml(pattern, m, store_to_file_rev);

        lens_rev.close();
        col_ids_rev.close();
    }

    void _query_pml(const char *pattern, const size_t m, std::function<void(const ulint, const ulint, const ulint)> store_stats)
    {
        std::vector<ulint> pml_lengths(m);
        std::vector<ulint> col_stats(m);

        ulint pos = n - 1;
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
        ulint thr = n;

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
#endif /* end of include guard: _TEST_BWT_HH */