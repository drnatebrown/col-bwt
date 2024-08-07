/* Lf_table - Uncompressed version of OptBWTR (LF table)
    Copyright (C) 2021 Nathaniel Brown
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
   \file LF_table.hpp
   \brief LF_table.hpp Uncompressed version of OptBWTR (LF table)
   \author Nathaniel Brown
   \date 19/11/2021
*/

#ifndef _LF_TABLE_HH
#define _LF_TABLE_HH

#include <common.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

using namespace std;

// Row of the LF table
class LF_row
{
public:
    char character;
    ulint idx : BWT_BITS;
    ulint interval : BWT_BITS;
    ulint offset : BWT_BITS;

    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name ="")
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        out.write((char *)&character, sizeof(character));
        written_bytes += sizeof(character);

        size_t temp = idx;
        out.write((char *)&temp, BWT_BYTES);
        written_bytes += BWT_BYTES;

        temp = interval;
        out.write((char *)&temp, BWT_BYTES);
        written_bytes += BWT_BYTES;

        temp = offset;
        out.write((char *)&temp, BWT_BYTES);
        written_bytes += BWT_BYTES;

        return written_bytes;
    }

    void load(std::istream &in)
    {
        in.read((char *)&character, sizeof(character));

        size_t temp;
        in.read((char *)&temp, BWT_BYTES);
        idx = temp;

        in.read((char *)&temp, BWT_BYTES);
        interval = temp;

        in.read((char *)&temp, BWT_BYTES);
        offset = temp;
    }
};

template <typename T = LF_row>
class LF_table
{
public:
    LF_table() {}

    LF_table(std::ifstream &heads, std::ifstream &lengths)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);
        
        LF_runs = vector<LF_row>();
        vector<vector<size_t>> L_block_indices = vector<vector<size_t>>(ALPHABET_SIZE);
        
        char c;
        ulint i = 0;
        r = 0;
        n = 0;
        while ((c = heads.get()) != EOF)
        {
            size_t length = 0;
            lengths.read((char *)&length, 5);
            if (c <= TERMINATOR) c = TERMINATOR;

            LF_runs.push_back({c, n, 0, 0});
            L_block_indices[c].push_back(i++);
            n+=length;
        }
        r = LF_runs.size();

        compute_table(L_block_indices);

        #ifdef PRINT_STATS
        cout << "Text runs: " << runs() << std::endl;
        cout << "Text length: " << size() << std::endl;
        #endif
    }

    LF_table(std::ifstream &bwt)
    {
        bwt.clear();
        bwt.seekg(0);
        
        LF_runs = vector<LF_row>();
        vector<vector<size_t>> L_block_indices = vector<vector<size_t>>(ALPHABET_SIZE);
        
        char last_c;
        char c;
        ulint i = 0;
        r = 0;
        n = 0;
        size_t length = 0;
        while ((c = bwt.get()) != EOF)
        {
            if (c <= TERMINATOR) c = TERMINATOR;
            
            if (i != 0 && c != last_c)
            {
                LF_runs.push_back({last_c, n, 0, 0});
                L_block_indices[last_c].push_back(i++); 
                n+=length;
                length = 0;
            }
            ++length;
            last_c = c;
        }
        // Step for final character
        LF_runs.push_back({last_c, length, 0, 0});
        L_block_indices[last_c].push_back(i++);  
        n+=length;

        r = LF_runs.size();

        compute_table(L_block_indices);

        #ifdef PRINT_STATS
        cout << "Text runs: " << runs() << std::endl;
        cout << "Text length: " << size() << std::endl;
        #endif
    }

    const LF_row get(size_t i)
    {
        assert(i < LF_runs.size());
        return LF_runs[i];
    }

    uchar get_char(ulint i)
    {
        return get(i).character;
    }

    // TODO specilization for template with lengths, not idx
    ulint get_length(ulint i)
    {
        return (i == r - 1) ? n : (get_idx(i + 1) - get_idx(i));
    }

    ulint get_idx(ulint i)
    {
        return get(i).idx;
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

        char c;
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

    std::string get_file_extension() const
    {
        return ".LF_table";
    }

    void bwt_stats()
    {
        cout << "Number of BWT equal-letter runs: r = " <<  r << std::endl;
        cout << "Length of complete BWT: n = " << n << std::endl;
        cout << "Rate n/r = " << double(n) / r <<  r << std::endl;
        cout << "log2(r) = " << log2(double(r)) <<  r << std::endl;
        cout << "log2(n/r) = " << log2(double(n) / r) <<  r << std::endl;
    }

    void mem_stats()
    {
        cout << "Memory consumption (bytes)." << std::endl;
        cout << "   LF Table: " << r + 4*r*BWT_BYTES << std::endl;
        cout << "           Chars: " << r << std::endl;
        cout << "            Ints: " << r*BWT_BYTES << std::endl;
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

        size_t size = LF_runs.size();
        out.write((char *)&size, sizeof(size));
        written_bytes += sizeof(size);

        for(size_t i = 0; i < size; ++i)
        {
            written_bytes += LF_runs[i].serialize(out, v, "LF_run_" + std::to_string(i));
        }

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

        in.read((char *)&size, sizeof(size));
        LF_runs = std::vector<T>(size);
        for(size_t i = 0; i < size; ++i)
        {
            LF_runs[i].load(in);
        }
    }

protected:
    ulint n; // Length of BWT
    ulint r; // Runs of BWT

    vector<T> LF_runs;

    void compute_table(vector<vector<ulint>> L_block_indices) {
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
};

#endif /* end of include guard: _LF_TABLE_HH */