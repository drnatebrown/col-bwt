#include <common.hpp>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <malloc_count.h>
#include <sys/stat.h>

int main(int argc, char *const argv[]) {
    Args args;
    parseArgs(argc, argv, args);

    message("Creating BWT from RLBWT");
    std::string head_filename = args.filename + ".bwt.heads";
    std::string len_filename = args.filename + ".bwt.len";
    std::string output_filename = args.filename + ".bwt";

    std::ifstream head_file(head_filename, std::ios::binary);
    std::ifstream len_file(len_filename, std::ios::binary);
    std::ofstream output_file(output_filename, std::ios::binary);

    status("Writing Characters to BWT");
    char bwt_char;
    uint64_t bwt_len;
    while (head_file.get(bwt_char) && len_file.read(reinterpret_cast<char*>(&bwt_len), RW_BYTES)) {
        std::string output_string(bwt_len, bwt_char);
        output_file.write(output_string.c_str(), output_string.size());
    }
    status();

    head_file.close();
    len_file.close();
    output_file.close();
    return 0;
}