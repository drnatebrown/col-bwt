#include <common.hpp>
#include <col_bwt.hpp>

#include <sdsl/io.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>

int main(int argc, char *const argv[])
{
    Timer timer;
    Args args;
    parseArgs(argc, argv, args);

    message("Building Col BWT supporting Co-linearity Statistics: ");
    timer.start();

    std::string col_ids_fname = args.filename + ".col_ids";
    std::ifstream ifs_col_ids(col_ids_fname);
    std::string col_runs_fname = args.filename + ".col_runs";
    std::ifstream ifs_col_runs(col_runs_fname);
    std::string thr_fname = args.filename + ".thr_pos";
    std::ifstream ifs_thr(thr_fname);

    sdsl::sd_vector<> col_runs;
    col_runs.load(ifs_col_runs);

    std::string bwt_fname = args.filename + ".bwt";

    std::string bwt_heads_fname = bwt_fname + ".heads";
    std::ifstream ifs_heads(bwt_heads_fname);
    std::string bwt_len_fname = bwt_fname + ".len";
    std::ifstream ifs_len(bwt_len_fname);

    ifs_heads.seekg(0);
    ifs_len.seekg(0);
    ifs_col_ids.seekg(0);
    // col_bwt col_bwt(ifs_heads, ifs_len, ifs_col_ids, col_runs);
    col_pml tbl(ifs_heads, ifs_len, ifs_col_ids, ifs_thr, col_runs);
    timer.end();

    tbl.bwt_stats();
    tbl.mem_stats();

    submessage("Construction Complete");
    timer.startTime();

    message("Serializing");
    timer.mid();

    std::string tbl_outfile = args.filename + tbl.get_file_extension();
    std::ofstream out_fp(tbl_outfile);
    tbl.serialize(out_fp);
    out_fp.close();

    timer.end();
    submessage("Serializing Complete");
    timer.midTime();

    message("Done", false);
    mem_peak();
    timer.startTime();

    return 0;
}