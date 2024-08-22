/* rle_string_w - RLBWT based on r-index rle_string
    Copyright (C) 2021 Massimiliano Rossi

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
   \file rle_string_w.hpp
   \brief rle_string_w.hpp RLBWT based on r-index rle_string.
   \author Massimiliano Rossi
   \date 04/11/2021
*/


#ifndef _RLE_STRING_W_HH
#define _RLE_STRING_W_HH

#include <iostream>
#include <rle_string.hpp>

using namespace ri;
// Lifting from the r-index source code of the RLBWT
template <class rle_string_t = ri::rle_string_sd>
class r_index {
public :

  static const uchar TERMINATOR = 1;
  rle_string_t bwt;
  vector<ulint> F;
  ulint terminator_position = 0;

  typedef size_t size_type;

  r_index(){}

  r_index(std::string filename)
  {
        std::string bwt_fname = filename + ".bwt";

        std::string bwt_heads_fname = bwt_fname + ".heads";
        std::ifstream ifs_heads(bwt_heads_fname);
        assert(ifs_heads.is_open());
        std::string bwt_len_fname = bwt_fname + ".len";
        std::ifstream ifs_len(bwt_len_fname);
        assert(ifs_len.is_open());
        bwt = rle_string_t(ifs_heads, ifs_len);

        ifs_heads.seekg(0);
        ifs_len.seekg(0);
        build_F_(ifs_heads, ifs_len);
  }

  /*
  * \param i position in the BWT
  * \param c character
  * \return lexicographic rank of cw in bwt
  */
  ulint LF(ulint i, uchar c)
  {
      //number of c before the interval
      ulint c_before = bwt.rank(i, c);
      // number of c inside the interval rn
      ulint l = F[c] + c_before;
      return l;
  }

  /*
    * \param r inclusive range of a string w
    * \param c character
    * \return inclusive range of cw
    */
  range_t LF(range_t rn, uchar c){
      //if character does not appear in the text, return empty pair
      if((c==255 and F[c]==bwt.size()) || F[c]>=F[c+1])
          return {1,0};
      //number of c before the interval
      ulint c_before = bwt.rank(rn.first,c);
      // number of c inside the interval rn
      // if rn.1 == rn.2, then this redudant rank query _should_ be avoided,
      // but this case won't be encountered too often if n/r is high
      ulint c_inside = bwt.rank(rn.second+1,c) - c_before;
      //if there are no c in the interval, return empty range
      if(c_inside==0) return {1,0};
      ulint l = F[c] + c_before;
      return {l,l+c_inside-1};
  }

    //forward navigation of the BWT
    uint64_t FL(uint64_t  i) {
        //i-th character in first BWT column
        auto c = f_at(i);
        //this c is the j-th (counting from 0)
        uint64_t j = i - F[c];
        return bwt.select(j,uint8_t(c));
    }

    //forward navigation of the BWT, where for efficiency we give c=F[i] as input
    // uint64_t FL(uint64_t  i, uint8_t c) {
    //     //i-th character in first BWT column
    //     assert(c == f_at(i));
    //     //this c is the j-th (counting from 0)
    //     uint64_t j = i - F[c];
    //     return bwt.select(j,uint8_t(c));
    // }

    //forward navigation of the BWT using interval/offset representation
    uint64_t FL(uint64_t k, uint64_t d) {
        return FL(run_start(k) + d);
    }

  /* serialize the structure to the ostream
  * \param out     the ostream
  */
  size_t serialize(std::ostream& out){
      ulint w_bytes = 0;
      assert(F.size()>0);
      assert(bwt.size()>0);
      out.write((char*)&terminator_position,sizeof(terminator_position));
      out.write((char*)F.data(),256*sizeof(ulint));
      w_bytes += sizeof(terminator_position) + 256*sizeof(ulint);
      w_bytes += bwt.serialize(out);
      return w_bytes;
  }

  /* load the structure from the istream
  * \param in the istream
  */
  void load(std::istream& in) {
      in.read((char*)&terminator_position,sizeof(terminator_position));
      F = vector<ulint>(256);
      in.read((char*)F.data(),256*sizeof(ulint));
      bwt.load(in);
  }

  std::string get_file_extension() const
  {
      return ".rle";
  }

  size_t size()
  {
      return bwt.size();
  }

   size_t runs()
  {
      return bwt.number_of_runs();
  }

  size_t run_start(size_t j) {
    return bwt.run_range(j).first;
  }

  ulint get_run(ulint i)
  {
      return bwt.run_of_position(i);
  }

  ulint get_length(ulint k)
  {
      return bwt.run_at(k);
  }

  uchar get_char(size_t i)
  {
    return f_at(i);
  }

protected:
  vector<ulint> build_F_(std::ifstream &heads, std::ifstream &lengths)
  {
      heads.clear();
      heads.seekg(0);
      lengths.clear();
      lengths.seekg(0);

      F = vector<ulint>(256, 0);
      int c;
      ulint i = 0;
      while ((c = heads.get()) != EOF)
      {
          size_t length = 0;
          lengths.read((char *)&length, 5);
          if (c > TERMINATOR)
              F[c] += length;
          else
          {
              F[TERMINATOR] += length;
              terminator_position = i;
          }
          i++;
      }
      for (ulint i = 255; i > 0; --i)
          F[i] = F[i - 1];
      F[0] = 0;
      for (ulint i = 1; i < 256; ++i)
          F[i] += F[i - 1];
      return F;
  }

  uint8_t f_at(uint64_t i) const {
        uint64_t c = (std::upper_bound(F.begin(), F.end(),i) - F.begin()) - 1;
        return uint8_t(c);
  }

};

#endif /* end of include guard: _RLE_STRING_W_HH */
