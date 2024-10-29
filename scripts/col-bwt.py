#!/usr/bin/env python3

import argparse
import os
import subprocess

MUMEMTO_DIR = "_deps/mumemto-build/"
MOVI_DIR = "_deps/movi-build/"
COL_BWT_DIR = "./src/"

mumemto_prg = os.path.join(MUMEMTO_DIR, "mumemto")
rlbwt_to_bwt_prg = os.path.join(COL_BWT_DIR, "rlbwt_to_bwt")
col_split_build = os.path.join(COL_BWT_DIR, "build_FL")
col_split_prg = os.path.join(COL_BWT_DIR, "col_split")
movi_split = os.path.join(MOVI_DIR, "movi-split")

ascii_art = """                                              
         /\                 /\                               
        /**\           /\  /**\                              
       /****\   /\    /**\/****\                             
      /      \ /**\  /   /      \                            
     /  /\    /    \/   /        \ _       _             _   
    /  /  \  /      \  /  ___ ___ | |     | |____      _| |_ 
   /  /    \/ /\     \/  / __/ _ \| |_____| '_ \ \ /\ / / __|
  /  /      \/  \/\   \ | (_| (_) | |_____| |_) \ V  V /| |_ 
_/__/_______/___/__\___\ \___\___/|_|     |_.__/ \_/\_/  \__|
"""

def file_exists(filename):
    return os.path.exists(filename)

def exec_command(command):
    print(f"\033[1;33m>> {' '.join(command)}\033[0m")
    result = subprocess.run(["/usr/bin/time", "--format=Wall Time: %e\nMax Memory: %M"] + command)
    print(f"\033[1;33m<<\033[0m")
    return result.returncode

def big_log(message):
    print(f"\033[1;31m{message}\033[0m")

def log(*args):
    if len(args) == 2:
        print(f"\033[1;32m{args[0]} \033[1;35m{args[1]}\033[0m")
    else:
        print(f"\033[1;35m{args[0]}\033[0m")

def build(args):
    filename = f"{args.output}.fa"
    mum_file = f"{filename}.col_mums"

    head_file = f"{filename}.bwt.heads"
    len_file = f"{filename}.bwt.len"
    bwt_file = f"{filename}.bwt"
    thr_file = f"{filename}.thr_pos"
    parse_file = f"{filename}.parse"
    dict_file = f"{filename}.dict"

    big_log(f"[[ Processing {filename} ]]")

    N = int(subprocess.check_output(["wc", "-l", args.input]).split()[0])
    big_log(f"Number of sequences: {N}")

    flags = []
    if not file_exists(head_file) or not file_exists(len_file) or args.force:
        flags.append("-R")
    if not file_exists(thr_file) or args.force:
        flags.append("-T")
    if file_exists(mum_file) and len(flags) == 0 and not args.force:
        log("[mumemto]", "Multi-MUMs already found, skipping...")
    else:
        if file_exists(parse_file) and file_exists(dict_file):
            flags.append("-s")
        log("[mumemto]", f"Running mumemto to generate multi-MUMs")
        status = exec_command([mumemto_prg, "mum", "-i", args.input, "-o", args.output, "-r", "-K", "-l", f"{args.min_mum}"] + flags)
        if status != 0:
            log("[mumemto]", "Error running mumemto, exiting...")
            subprocess.run(["rm", "-f", mum_file, head_file, len_file,thr_file, parse_file, dict_file])
            sys.exit(status)
    
    flags = []
    if args.verbose:
        flags.append("-v")

    FL_table = f"{filename}.FL_table"
    if file_exists(FL_table) and not args.force:
        log("[col-bwt]", "FL move structure already exists, skipping...")
    else:
        log("[col-bwt]", f"Building FL move structure")
        status = exec_command([col_split_build, filename] + flags)
        if status != 0:
            log("[col-bwt]", "Error building FL table, exiting...")
            subprocess.run(["rm", "-f", FL_table])
            sys.exit(status)
    
    log("[col-bwt]",f"Finding multi-MUM sub-runs with mode {args.mode} and split rate {args.sub_sample}")
    status = exec_command([col_split_prg, filename, "-m", args.mode, "-s", f"{args.sub_sample}"] + flags)
    if status != 0:
        log("[col-bwt]", "Error splitting RLBWT, exiting...")
        sys.exit(status)

    if not file_exists(f"{filename}.bwt") and not args.force:
        log("[col-bwt]", "Generating BWT for Movi")
        status = exec_command([rlbwt_to_bwt_prg, filename] + flags)
        if status != 0:
            log("[col-bwt]", "Error generating BWT, exiting...")
            subprocess.run(["rm", "-f", bwt_file])
            sys.exit(status)

    log("[movi]", "Building col-bwt")
    status = exec_command([movi_split, "build", "-f", filename, "-i", args.output] + flags)
    if status != 0:
        log("[movi]", "Error building col-bwt, exiting...")
        sys.exit(status)

    big_log(f"[[ Finished processing {filename} ]]")

def query(args):
    big_log(f"Generating PMLs and chain statistics for {args.pattern}")
    log("[movi]", "Querying col-bwt")
    status = exec_command([movi_split, "query", "-i", args.index, "-r", args.pattern])
    if status != 0:
        log("[col-bwt]", "Error running Movi, exiting...")
        sys.exit(status)


def main():
    parser = argparse.ArgumentParser(description="Full-text index for pangenomes using chain statistics")
    subparsers = parser.add_subparsers(dest="command")

    # Build command
    build_parser = subparsers.add_parser("build", help="Finds multi-MUMs and uses to build col-bwt")
    build_parser.add_argument("-i", "--input", required=True, help="path to a file-list of genomes to use (overrides positional args)", type=str)
    build_parser.add_argument("-o", "--output", required=True, help="output prefix path", type=str)
    build_parser.add_argument("-m", "--mode", 
                             help="splitting mode (options: tunnels, all)", type=str, default="tunnels")
    build_parser.add_argument("-s", "--sub-sample", 
                             help="sub-sample rate for splitting", type=int, default=10)
    build_parser.add_argument("-l", "--min-mum", help="minimum multi-MUM length", type=int, default=20)
    build_parser.add_argument("-v", "--verbose", action="store_true", help="verbose output")
    build_parser.add_argument("--force", action="store_true", help="forces all build steps to run")
    # build_parser.add_argument("--keep", action="store_true", help="keep temporary files")
    # add keep
    # add n args (like mumemto)
    # build_parser.add_argument("-n", "--n-args", nargs='+', help="Additional arguments for mumemto")

    # Query command
    query_parser = subparsers.add_parser("query", help="Computes PMLs and chain statistics for pattern")
    query_parser.add_argument("index", type=str, help="output prefix path of build command")
    query_parser.add_argument("-p", "--pattern", required=True, help="pattern fasta file")
    # query_parser.add_argument("--movi", help="access the Movi exe directly, ignoring col-bwt flags, action="store_true")

    args = parser.parse_args()

    if args.command == "build":
        build(args)
    elif args.command == "query":
        query(args)
    else:
        print(ascii_art)
        parser.print_help()

if __name__ == "__main__":
    main()