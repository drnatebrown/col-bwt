/*  common
    Copyright (C) 2020 Massimiliano Rossi
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
   \file common.hpp
   \brief common.hpp contains common features.
   \author Massimiliano Rossi
   \date 12/03/2020
*/

#ifndef _COMMON_HH
#define _COMMON_HH

#include <cstdlib>
#include <cstdio>
#include <ctime>

#include <sys/time.h>

#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <vector>      // std::vector
#include <assert.h>

#include <sdsl/int_vector.hpp>

#define ALPHABET_SIZE 256
#define RW_BYTES 5
#define BWT_BYTES 5
#define PRINT_STATS 1

static const uint8_t TERMINATOR = 1;
typedef unsigned long int ulint;
typedef unsigned char uchar;

//*********************** Argument options ***************************************
// struct containing command line parameters and other globals
struct Args
{
  std::string filename = "";
  bool rle   = true; // read in RLBWT
  int N = 0; // size of collection
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
  int c;
  extern char *optarg;
  extern int optind;

  std::string sarg;
  while ((c = getopt(argc, argv, "rd:N:")) != -1)
  {
    switch (c)
    {
    case 'r':
      arg.rle = true;
      break;
    case 'N':
      sarg.assign(optarg);
      arg.N = std::stoi(sarg);
      break;
    case '?':
      std::cout << "ERROR: Unknown option.\n";
      break;
    }
  }
  // the only input parameter is the file name
  if (argc == optind + 1)
  {
    arg.filename.assign(argv[optind]);
  }
  else
  {
    std::cout << "ERROR: Invalid number of arguments\n";
  }
}

//********** end argument options ********************

// Convert boolean vector to specified bit vector
template<class B>
B bool_to_bit_vec(std::vector<bool> &b)
{
  if(b.size()==0) return B();

  sdsl::bit_vector bv(b.size());

  for(size_t i = 0; i < b.size(); ++i)
    bv[i] = b[i];

  return B(bv);
}

uint8_t bitsize(uint64_t x){

	if(x==0) return 1;
	return 64 - __builtin_clzll(x);

}

#endif /* end of include guard: _COMMON_HH */