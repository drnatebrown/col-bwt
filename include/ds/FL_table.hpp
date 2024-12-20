/* FL_table - Table supporting first-to-last mapping of BWT
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
   \file FL_table.hpp
   \brief FL_table.hpp Table supporting first-to-last mapping of BWT
   \author Nathaniel Brown
   \author Massimiliano Rossi
   \date 07/08/2021
*/

#ifndef _FL_TABLE_HH
#define _FL_TABLE_HH

#include <common.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

using namespace std;

//TODO move to using reduced alphabet
class FL_table
{
public:
    // Row of the FL table
    struct FL_row
    {
        char character;
        ulint idx : BWT_BITS;
        ulint interval : RUN_BITS;
        ulint offset : LEN_BITS;

        size_t serialize(std::ostream &out)
        {
            size_t written_bytes = 0;

            out.write((char *)&character, sizeof(character));
            written_bytes += sizeof(character);

            size_t temp = idx;
            out.write((char *)&temp, BWT_BYTES);
            written_bytes += BWT_BYTES;

            temp = interval;
            out.write((char *)&temp, RUN_BYTES);
            written_bytes += RUN_BYTES;

            temp = offset;
            out.write((char *)&temp, LEN_BYTES);
            written_bytes += LEN_BYTES;

            return written_bytes;
        }

        void load(std::istream &in)
        {
            in.read((char *)&character, sizeof(character));

            size_t temp;
            in.read((char *)&temp, BWT_BYTES);
            idx = temp;

            in.read((char *)&temp, RUN_BYTES);
            interval = temp;

            in.read((char *)&temp, LEN_BYTES);
            offset = temp;
        }
    } __attribute__((packed));


    FL_table() {}

    FL_table(std::ifstream &heads, std::ifstream &lengths)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);
        
        vector<char> L_chars = vector<char>();
        vector<ulint> L_lens = vector<ulint>();
        vector<vector<ulint>> L_block_indices = vector<vector<ulint>>(256);
        vector<vector<ulint>> char_runs = vector<vector<ulint>>(256); // Vector containing lengths for runs of certain character

        status("Constructing BWT table (FL)");
        char c = 0;
        n = 0;
        ulint i = 0;
        while ((c = heads.get()) != EOF)
        {
            if (c <= TERMINATOR) c = TERMINATOR;

            size_t length = 0;
            lengths.read((char *)&length, RW_BYTES);

            L_chars.push_back(c);
            L_lens.push_back(length);
            L_block_indices[c].push_back(i++);
            char_runs[c].push_back(length);

            n+=length;
        }
        r = L_chars.size();
        status();

        status("Computing values of BWT table");
        compute_table(L_chars, L_lens, L_block_indices, char_runs);
        status();

        status("Computing bitvector marking run-heads in L");
        compute_L_heads(L_lens);
        status();

        #ifdef PRINT_STATS
        stat("Text runs", runs());
        stat("Text length", size());
        #endif
    }

    FL_table(std::ifstream &bwt)
    {
        bwt.clear();
        bwt.seekg(0);
        
        vector<char> L_chars = vector<char>();
        vector<ulint> L_lens = vector<ulint>();
        vector<vector<ulint>> L_block_indices = vector<vector<ulint>>(256);
        vector<vector<ulint>> char_runs = vector<vector<ulint>>(256); // Vector containing lengths for runs of certain character
        
        status("Constructing BWT table (FL)");
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
                L_chars.push_back(c);
                L_lens.push_back(length);
                L_block_indices[c].push_back(i++);
                char_runs[c].push_back(length);
                n+=length;
                length = 0;
            }
            ++length;
            last_c = c;
        }
        // Step for final character
        L_chars.push_back(c);
        L_lens.push_back(length);
        L_block_indices[c].push_back(i++);
        char_runs[c].push_back(length);
        n+=length;

        r = L_chars.size();
        status();

        status("Computing values of BWT table");
        compute_table(L_chars, L_lens, L_block_indices, char_runs);
        status();

        status("Computing bitvector marking run-heads in L");
        compute_L_heads(L_lens);
        status();

        #ifdef PRINT_STATS
        stat("Text runs", runs());
        stat("Text length", size());
        #endif
    }

    const FL_row get(size_t i)
    {
        assert(i < FL_runs.size());
        return FL_runs[i];
    }

    ulint size()
    {
        return n;
    }

    ulint runs()
    {
        return r;
    }

    void decompress(std::string outfile) 
    {
        std::ofstream out(outfile);

        // Forward step to first character
        std::pair<ulint, ulint> pos = FL(0, 0);
        pos = FL(pos.first, pos.second);

        char c;
        while((c = get_char(pos.first)) > TERMINATOR) 
        {
            out << c;
            pos = FL(pos.first, pos.second);
        }
    }

    /*
     * \param Run position (RLE intervals)
     * \param Current character offset in block
     * \return block position and offset of preceding character
     */
    std::pair<ulint, ulint> FL(ulint run, ulint offset)
    {
        ulint next_interval = FL_runs[run].interval;
	    ulint next_offset = FL_runs[run].offset + offset;

	    while (next_offset >= get_length(next_interval)) 
        {
            next_offset -= get_length(next_interval++);
        }

	    return std::make_pair(next_interval, next_offset);
    }

    uchar get_char(ulint i)
    {
        return get(i).character;
    }

    ulint get_length(ulint i)
    {
        return (i == r - 1) ? (n - get_idx(i)) : (get_idx(i + 1) - get_idx(i));
    }

    ulint get_idx(ulint i)
    {
        return get(i).idx;
    }

    std::string get_file_extension() const
    {
        return ".FL_table";
    }

    ulint get_L_pos(ulint i)
    {
        return L_heads_select(i + 1);
    }

    ulint get_L_length(ulint i)
    {
        return (i == r - 1) ? (n - get_L_pos(i)) : (get_L_pos(i + 1) - get_L_pos(i));
    }

    bool is_L_head(ulint i)
    {
        return L_heads[i];
    }

    sdsl::sd_vector<>& get_L_heads()
    {
        return L_heads;
    }

    void bwt_stats()
    {
        log("Number of BWT equal-letter runs: r = ", r);
        log("Length of complete BWT: n = ", n);
        log("Rate n/r = ", double(n) / r);
        log("log2(r) = ", log2(double(r)));
        log("log2(n/r) = ", log2(double(n) / r));
    }

    void mem_stats()
    {
        sdsl::nullstream ns;
        size_t L_runs_size = L_heads.serialize(ns);
        log("Memory (bytes):");
        log("   FL Table: ", sizeof(FL_row)*r);
        log("    L Heads: ", L_runs_size);
        log("       Ints: ", sizeof(n) + sizeof(r));
        log("      Total: ", sizeof(FL_row)*r + L_runs_size + sizeof(n) + sizeof(r));
    }

    /* serialize to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name ="")
    {
        size_t written_bytes = 0;

        out.write((char *)&n, sizeof(n));
        written_bytes += sizeof(n);

        out.write((char *)&r, sizeof(r));
        written_bytes += sizeof(r);

        written_bytes += L_heads.serialize(out);

        written_bytes += write_vec(FL_runs, out);

        return written_bytes;
    }

    /* load from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        in.read((char *)&n, sizeof(n));
        in.read((char *)&r, sizeof(r));

        L_heads.load(in);
        L_heads_select = sd_select(&L_heads);

        FL_runs = std::vector<FL_row>(r);
        read_vec(FL_runs, in);
    }
    
protected:
    ulint n; // Length of BWT
    ulint r; // Runs of BWT

    vector<FL_row> FL_runs;
    sdsl::sd_vector<> L_heads;
    sd_select L_heads_select;

    void compute_table(vector<char> L_chars, vector<ulint> L_lens, vector<vector<ulint>> L_block_indices, vector<vector<ulint>> char_runs) {
        FL_runs = vector<FL_row>(r);
        size_t i = 0;
        size_t curr_idx = 0;
        for (size_t c = 0; c < 256; ++c)
        {
            for (size_t j = 0; j < char_runs[c].size(); ++j) {
                size_t length = char_runs[c][j];
                FL_runs[i].character = (unsigned char) c;
                FL_runs[i].idx = curr_idx;
                ++i;
                curr_idx += length;
            }
        }

        ulint k = 0; // current row to be filled
        for(size_t i = 0; i < L_block_indices.size(); ++i) 
        {
            ulint F_curr = 0; // current position when scanning F
            ulint F_seen = 0; // characters seen before position in F
            ulint L_curr = 0; // current position when scanning L
            ulint L_seen = 0; // characters seen before position in L
            for(size_t j = 0; j < L_block_indices[i].size(); ++j) 
            {
                while (L_curr < L_block_indices[i][j]) {
                    L_seen += L_lens[L_curr++];
                }
                while (F_seen + get_length(F_curr) <= L_seen) {
                    F_seen += get_length(F_curr++);
                }

                FL_runs[k].interval = F_curr;
                FL_runs[k].offset = L_seen - F_seen;
                ++k;
            }
        }
    }

    void compute_L_heads(vector<ulint> &L_lens) {
        assert(L_lens.size() == r);
        sdsl::sd_vector_builder L_heads_builder(n, L_lens.size());
        size_t idx = 0;
        for (ulint i = 0; i < L_lens.size(); ++i) {
            L_heads_builder.set(idx);
            idx += L_lens[i];
        }
        L_heads = sdsl::sd_vector<>(L_heads_builder);
        L_heads_select = sd_select(&L_heads);
    }
};

#endif /* end of include guard: _FL_TABLE_HH */