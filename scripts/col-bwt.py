#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess

MUMEMTO_DIR = "_deps/mumemto-build/"
MOVI_DIR = "_deps/movi-build/"
COL_BWT_DIR = "./src/"

mumemto_prg = os.path.join(MUMEMTO_DIR, "mumemto")
rlbwt_to_bwt_prg = os.path.join(COL_BWT_DIR, "rlbwt_to_bwt")
FL_build_prg = os.path.join(COL_BWT_DIR, "build_FL")
col_split_prg = os.path.join(COL_BWT_DIR, "col_split")
movi_split = os.path.join(MOVI_DIR, "movi-split")

# files generated during runtime
mum_temp_exts = ["last", "occ", "parse_old", "parse", "dict"]
col_temp_exts = ["bwt.heads", "bwt.len", "col_ids", "col_runs"]
clean_exts = ["bwt", "thr_pos", "FL_table", "col_mums"] + mum_temp_exts + col_temp_exts

# use different file base
prefix_exts = ["fa", "lengths"]

ascii_art = """                                              
         /\                 /\                          /\     
        /**\           /\  /**\                        /**\    
       /****\   /\    /**\/****\                    __/****\     
      /      \ /**\__/   /      \                  /        \ 
     /  /\    /    \    /        \ _       _      /      _   \__
    /  /  \  /      \  /  ___ ___ | |     | |____/     _| |     \ 
   /  /    \/ /\     \/  / __/ _ \| |_____| '_ \ \ /\ / / __|    \ 
  /  /      \/  \/\   \ | (_| (_) | |_____| |_) \ V  V /| |_      \ 
_/__/_______/___/__\___\_\___\___/|_|     |_.__/ \_/\_/  \__|______\_
"""

class OnError:
    def __init__(self, program="ERROR", message="Error running command", artifacts=[]):
        self.program = program
        self.message = message
        self.artifacts = artifacts

def file_exists(filename):
    return os.path.exists(filename)

def gather_files(filename, extensions):
    return [f"{filename}.{ext}" for ext in extensions]

def remove_files(files, force = True):
    if force:
        subprocess.run(["rm", "-f"] + files)
    else:
        subprocess.run(["rm"] + files)

def get_N(filename):
    try:
        result = subprocess.check_output(f"awk '{{if($2 > max) max=$2}} END {{print max}}' {filename}", shell=True)
        return int(result.strip())
    except subprocess.CalledProcessError as e:
        log_error(f"Input file not structured correctly")
        sys.exit(1)

def exec_command(command, on_error=OnError()):
    print(f"\033[1;33m>> {' '.join(command)}\033[0m")
    result = subprocess.run(["/usr/bin/time", "--format=Wall Time: %e\nMax Memory: %M"] + command)
    print(f"\033[1;33m<<\033[0m")
    check_status(result.returncode, on_error)

def check_status(status, on_error):
    if status != 0:
        rm_message = "removing artifacts and " if on_error.artifacts else ""
        log_error(f"[{on_error.program}]", f"{on_error.message}, {rm_message}exiting...")
        if on_error.artifacts:
            log(f"Artifacts: {' '.join(on_error.artifacts)}")
            remove_files(on_error.artifacts)
        sys.exit(status)

def big_log(message):
    print(f"\033[1;32m{message}\033[0m")

def log(*args):
    if len(args) == 2:
        print(f"\033[1;32m{args[0]} \033[1;35m{args[1]}\033[0m")
    else:
        print(f"\033[1;35m{args[0]}\033[0m")

def log_error(*args):
    if len(args) == 2:
        print(f"\033[1;31m{args[0]} \033[1;35m{args[1]}\033[0m")
    else:
        print(f"\033[1;31m[ERROR] \033[1;35m{args[0]}\033[0m")

def build(args):
    filename = f"{args.output}.fa"
    mum_file = f"{filename}.col_mums"

    head_file = f"{filename}.bwt.heads"
    len_file = f"{filename}.bwt.len"
    thr_file = f"{filename}.thr_pos"
    parse_file = f"{filename}.parse"
    dict_file = f"{filename}.dict"

    bwt_file = f"{filename}.bwt"
    FL_table = f"{filename}.FL_table"

    big_log(f"[[ Processing {filename} ]]")

    if (args.input):
        if (not file_exists(args.input)):
            log_error(f"Input file {args.input} does not exist, exiting...")
            sys.exit(1)
        # N = int(subprocess.check_output(["wc", "-l", args.input]).split()[0])
        N = get_N(args.input)

    else:
        N = len(args.fastas)
    big_log(f"Number of sequences: {N}\n")

    flags = []
    files_created = [mum_file]
    need_rlbwt = not file_exists(FL_table) or not file_exists(f"{filename}.bwt")
    if need_rlbwt and (not file_exists(head_file) or not file_exists(len_file)) or args.force:
        flags.append("-R")
        files_created.append(head_file)
        files_created.append(len_file)
    if not file_exists(thr_file) or args.force:
        flags.append("-T")
        files_created.append(thr_file)
    if file_exists(mum_file) and len(flags) == 0 and not args.force:
        log("[mumemto]", "Multi-MUMs already found, skipping...")
    else:
        files_created.append(gather_files(args.output, prefix_exts))
        if file_exists(parse_file) and file_exists(dict_file):
            flags.append("-s")
        else:
            files_created.append(gather_files(filename, mum_temp_exts))
        if args.rev_comp:
            flags.append("-r")
        log("[mumemto]", f"Running mumemto to generate multi-MUMs")
        on_error = OnError("mumemto", "Error running mumemto", files_created)
        if args.input:
            exec_command([mumemto_prg, "mum", "-K", "-i", args.input, "-o", args.output, "-l", f"{args.min_mum}"] + flags, on_error)
        else:
            exec_command([mumemto_prg, "mum", "-K", "-o", args.output, "-l", f"{args.min_mum}"] + flags + args.fastas, on_error)
        if not args.keep_pfp:
            log(f"[mumemto]", "Removing temporary mumemto files")
            remove_files(gather_files(filename, mum_temp_exts) + gather_files(args.output, prefix_exts))

    flags = []
    files_created = []
    if args.verbose:
        flags.append("-v")

    files_created = [bwt_file]
    if not file_exists(f"{filename}.bwt") and not args.force:
        log("[col-bwt]", "Generating BWT for Movi")
        on_error = OnError("col-bwt", "Error generating BWT", files_created)
        exec_command([rlbwt_to_bwt_prg, filename] + flags, on_error)

    files_created = [FL_table]
    if file_exists(FL_table) and not args.force:
        log("[col-bwt]", "FL move structure already exists, skipping...")
    else:
        log("[col-bwt]", f"Building FL move structure")
        on_error = OnError("col-bwt", "Error building FL table", files_created)
        exec_command([FL_build_prg, filename] + flags, on_error)
        if not args.keep_col:
            log(f"[col-bwt]", "Removing temporary RLBWT files")
            remove_files(gather_files(filename, col_temp_exts))
    
    log("[col-bwt]",f"Finding multi-MUM sub-runs with mode {args.mode} and split rate {args.sub_sample}")
    on_error = OnError("col-bwt", "Error splitting RLBWT")
    exec_command([col_split_prg, filename, "-m", args.mode, "-s", f"{args.sub_sample}"] + flags)

    log("[movi]", "Building col-bwt")
    on_error = OnError("movi", "Error building col-bwt")
    status = exec_command([movi_split, "build", "-f", filename, "-i", args.output] + flags, on_error)

    if not args.keep_col:
        log(f"[col-bwt]", "Removing temporary col-bwt files")
        remove_files(gather_files(filename, col_temp_exts))

    big_log(f"\n[[ Finished processing {filename} ]]")
    log(f"Index output at {args.output}")

    if (args.clean):
        big_log(f"Removing intermediate files")
        remove_files(gather_files(filename, clean_exts) + gather_files(filename, mum_temp_exts) + gather_files(filename, col_temp_exts) + gather_files(args.output, prefix_exts))

def query(args):
    big_log(f"Generating PMLs and chain statistics for {args.pattern}")
    log("[movi]", "Querying col-bwt")
    files_created = [f"{args.pattern}.split.pml.bin", f"{args.pattern}.split.cid.bin"]
    on_error = OnError("movi", "Error querying col-bwt", files_created)
    status = exec_command([movi_split, "query", "-i", args.index, "-r", args.pattern], on_error)
    big_log(f"\n[[ Finished querying {args.pattern} ]]")
    log(f"Output at {args.pattern}.split.pml.bin and {args.pattern}.split.cid.bin")

def main():
    parser = argparse.ArgumentParser(description="Full-text index for pangenomes using chain statistics")
    subparsers = parser.add_subparsers(dest="command")

    # Build command
    build_parser = subparsers.add_parser("build", help="Finds multi-MUMs and uses to build col-bwt")
    
    # Either required
    build_parser.add_argument("fastas", help="fasta files to index", type=str, nargs="*")
    build_parser.add_argument("-i", "--input", help="Path to a file-list of genomes to use (overrides positional args)", type=str)

    build_parser.add_argument("-o", "--output", required=True, help="output prefix path", type=str)
    build_parser.add_argument("-r", "--rev_comp", action="store_true", help="include the reverse complement of the input sequences", default=False)
    build_parser.add_argument("-m", "--mode", 
                             help="splitting mode (options: tunnels, all)", type=str, default="tunnels")
    build_parser.add_argument("-s", "--sub-sample", 
                             help="sub-sample rate for splitting", type=int, default=10)
    build_parser.add_argument("-l", "--min-mum", help="minimum multi-MUM length", type=int, default=20)
    build_parser.add_argument("-v", "--verbose", action="store_true", help="verbose output")
    build_parser.add_argument("--force", action="store_true", help="forces all build steps to run")
    # build_parser.add_argument("--keep-pfp", action="store_true", help="keep temporary mumemto files")
    # build_parser.add_argument("--keep-col", action="store_true", help="keep temporary col-bwt files")
    build_parser.add_argument("--keep", action="store_true", help="keep all temporary files")
    build_parser.add_argument("--clean", action="store_true", help="clean up all intermediate files")

    # Query command
    query_parser = subparsers.add_parser("query", help="Computes PMLs and chain statistics for pattern")
    query_parser.add_argument("index", type=str, help="output prefix path of build command")
    query_parser.add_argument("-p", "--pattern", required=True, help="pattern fasta file")
    # query_parser.add_argument("--movi", help="access the Movi exe directly, ignoring col-bwt flags, action="store_true")

    args = parser.parse_args()

    if args.command == "build":
        args.keep_pfp = args.keep
        args.keep_col = args.keep
        if not args.fastas and not args.input:
            print("Error: Either positional arguments 'fastas' or the '-i/--input' option is required.")
            build_parser.print_help()
            sys.exit(1)
        build(args)
    elif args.command == "query":
        query(args)
    else:
        print(ascii_art)
        parser.print_help()

if __name__ == "__main__":
    main()