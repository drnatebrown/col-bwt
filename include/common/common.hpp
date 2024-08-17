/*  common
    Copyright (C) 2024 Massimiliano Rossi
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
   \author Nathaniel Brown
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

#include <utility>
#include <vector>      // std::vector
#include <assert.h>

#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>

#include <malloc_count.h>

/* CONFIGURABLE */
#define RW_BYTES 5
#define BWT_BYTES 5
#define ID_BITS 40
#define ASCII_SIZE 256
#define PRINT_STATS
#define MULTI_THREAD
#define DNA_ALPHABET

/* NOT CONFIGURABLE */
#ifdef DNA_ALPHABET
#define ALPHABET_SIZE 6
#define ALPHABET_BITS 3
#else
#define ALPHABET_SIZE ASCII_SIZE
#define ALPHABET_BITS 8
#endif

#define BWT_BITS (BWT_BYTES * 8)
#define ID_BYTES ((ID_BITS + 7) / 8)

bool verbose = false;

static const uint8_t TERMINATOR = 1;
typedef unsigned long int ulint;
typedef unsigned char uchar;

typedef typename sdsl::sd_vector<>::select_1_type sd_select;
typedef typename sdsl::sd_vector<>::rank_1_type sd_rank;

//*********************** Message options ***************************************

template <typename T>
std::string to_str(const T& value) {
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

/* ALWAYS DISPLAYED */
//**************************************************************
bool sub = false;
inline std::string indent() {
    return (sub) ? "\t" : "";
}

inline void error(const std::string& msg) {
    std::cerr << "[ERROR]: " << msg << std::endl;
}

inline void message(const std::string& msg, bool header = true) {
    sub = header;
    std::cout << "[INFO]: " << msg << std::endl;
}

inline void submessage(const std::string& msg) {
    sub = true;
    std::cout << indent() << "[INFO]: " << msg << std::endl;
}

inline void mem_peak() {
    std::cout << indent() << "[INFO]: Memory peak: " << malloc_count_peak() << std::endl;
}

#ifdef PRINT_STATS
template <typename T>
void stat(const std::string& label, const T& value) {
    std::cout << label << ": " << to_str(value) << std::endl;
}
#endif

/* VERBOSE ONLY */
//**************************************************************
template <typename T>
inline void log(const std::string& msg, const T& value) {
    if (verbose) {
        std::cout << indent() << "[LOG]: " << msg << to_str(value) << std::endl;
    }
}

inline void log(const std::string& msg) {
    if (verbose) {
        std::cout << indent() << "[LOG]: " << msg << std::endl;
    }
}

inline void status(const std::string& msg) {
    if (verbose) {
        std::cout << indent() << "[STATUS]: " << msg << "..." << std::flush;
    }
}

inline void status() {
    if (verbose) {
        std::cout << " DONE" << std::endl;
    }
}

class Timer {
public:
    void start() {
        t_start = std::chrono::high_resolution_clock::now();
    }

    void mid() {
        t_mid = std::chrono::high_resolution_clock::now();
    }

    void end() {
        t_end = std::chrono::high_resolution_clock::now();
    }

    // Method to print elapsed time from start to end
    void startTime() const {
        write_time(std::chrono::duration<double, std::ratio<1>>(t_end - t_start).count());
    }

    // Method to print elapsed time from mid to end
    void midTime() const {
        write_time(std::chrono::duration<double, std::ratio<1>>(t_end - t_mid).count());
    }

private:
    std::chrono::high_resolution_clock::time_point t_start;
    std::chrono::high_resolution_clock::time_point t_mid;
    std::chrono::high_resolution_clock::time_point t_end;

    void write_time(double duration) const {
        if (sub) {
              submessage("Elapsed time (s): " + to_str(duration));
          }
        else {
            message("Elapsed time (s): " + to_str(duration));
        }
    }
};

//********** end message options ********************

//*********************** Argument options ***************************************
// struct containing command line parameters and other globals
struct Args
{
    std::string filename = "";
    std::string pattern_filename = "";
    bool rle = true; // read in RLBWT
    bool verbose = false; // verbose output
    bool long_pattern = false; // write statistics direct to file
    int N = 0; // size of collection
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
    int c;
    extern char *optarg;
    extern int optind;

    std::string sarg;
    while ((c = getopt(argc, argv, "rvlN:p:")) != -1)
    {
        switch (c)
        {
        case 'r':
            arg.rle = true;
            break;
        case 'v':
            arg.verbose = true;
            verbose = arg.verbose;
            break;
        case 'l':
            arg.long_pattern = true;
            break;
        case 'N':
            sarg.assign(optarg);
            arg.N = std::stoi(sarg);
            break;
        case 'p':
            arg.pattern_filename.assign(optarg);
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
        error("Invalid number of arguments");
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

ulint bits_to_bytes(const ulint bits){
    return (bits + 7) / 8;
}

template <typename T>
size_t write_vec(const std::vector<T> &vec, std::ostream &out){
    const char* data = reinterpret_cast<const char*>(vec.data());
    size_t size = vec.size() * sizeof(T);
    out.write(data, size);
    return size;
}

template <typename T>
void read_vec(std::vector<T> &vec, std::istream &in){
    char* data = reinterpret_cast<char*>(vec.data());
    size_t size = vec.size() * sizeof(T);
    in.read(data, size);
}

#ifdef DNA_ALPHABET
constexpr std::array<uchar, ASCII_SIZE> initCharToBits() {
    std::array<uchar, ASCII_SIZE> table = {};
    table[TERMINATOR] = 0b000;
    table['A'] = 0b001;
    table['C'] = 0b010;
    table['G'] = 0b011;
    table['N'] = 0b100;
    table['T'] = 0b101;
    return table;
}
constexpr std::array<uchar, ASCII_SIZE> charToBits = initCharToBits();

constexpr uchar bitsToChar[ALPHABET_SIZE] = {
    TERMINATOR, // 000
    'A',  // 001
    'C',  // 010
    'G',  // 011
    'N',   // 100
    'T'  // 101
};
#endif

#endif /* end of include guard: _COMMON_HH */