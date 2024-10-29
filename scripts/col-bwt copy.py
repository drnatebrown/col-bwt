#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path

import os
import subprocess

MUMEMTO_DIR = "_deps/mumemto-build/src"
MOVI_DIR = "_deps/movi-build/src"
COL_BWT_DIR = "./src"

mumemto_prg = os.path.join(MUMEMTO_DIR, "build/mumemto")
col_split_build = os.path.join(COL_BWT_DIR, "build/src/build_FL")
col_split_prg = os.path.join(COL_BWT_DIR, "build/src/col_split")
rlbwt_to_bwt_prg = os.path.join(COL_BWT_DIR, "build/src/to_run_output")
movi_build_split = os.path.join(MOVI_DIR, "build/movi-split")

def Time(command):
    result = subprocess.run(["/usr/bin/time", "--format=Wall Time: %e\nMax Memory: %M"] + command, capture_output=True, text=True)
    print(result.stdout)
    print(result.stderr)

def exec_command(command):
    print(f"\033[1;33m>> {' '.join(command)}\033[0m")
    Time(command)
    print(f"\033[1;33m<<\033[0m")

def big_log(message):
    print(f"\033[1;31m{message}\033[0m")

def log(*args):
    if len(args) == 2:
        print(f"\033[1;32m{args[0]} \033[1;35m{args[1]}\033[0m")
    else:
        print(f"\033[1;35m{args[0]}\033[0m")

def prepare(args):
    big_log(f"[[ Processing {filename} ]]")

    N = int(subprocess.check_output(["wc", "-l", dataset_files]).split()[0])

    log("Number of sequences:", N)

    mum_file = f"{dataset_sep}.mums"
    if not os.path.exists(mum_file):
        log("[mumemto]", f"Building RLBWT, Thresholds; finding MUMs of length at least {os.getenv('MIN_MUM', '10')} on {filename}")
        pfp_flag = "-s" if os.path.exists(f"{dataset_sep}.parse") else ""
        bwt_flag = "-R" if not os.path.exists(f"{dataset_sep}.bwt.heads") else ""
        thr_flag = "-T" if not os.path.exists(f"{dataset_sep}.thr") else ""
        mumemto_filename = os.path.join(os.getenv("DATA_DIR", ""), filename)
        command = [
            mumemto_prg, "mum",
            "-i", dataset_files,
            "-o", mumemto_filename,
            "-r", "-K",
            "-l", os.getenv("MIN_MUM", "10"),
            pfp_flag, bwt_flag, thr_flag
        ]
        with open(log_file, "a") as log_f:
            subprocess.run(command, stdout=log_f, stderr=subprocess.STDOUT)
        os.rename(f"{mumemto_filename}.mums", mum_file)

    if not os.path.exists(f"{dataset_sep}.FL_table"):
        log("[bwt_frag]", f"Building FL table for {filename}")
        command = [col_split_build, dataset_sep]
        with open(log_file, "a") as log_f:
            subprocess.run(command, stdout=log_f, stderr=subprocess.STDOUT)

    if not os.path.exists(f"{dataset_sep}.bwt"):
        log("[bwt_frag]", f"Building raw BWT for {filename}")
        command = [rlbwt_to_bwt_prg, dataset_sep]
        with open(log_file, "a") as log_f:
            subprocess.run(command, stdout=log_f, stderr=subprocess.STDOUT)

    if not os.path.exists(f"{dataset_sep}.col_runs"):
        log("[bwt_frag]", f"Splitting RLBWT to create runs for MUMs of {filename} with tunnels and split rate 10")
        command = [col_split_prg, dataset_sep, "-N", str(N), "-m", "tunnels", "-v", "-s", "10"]
        with open(log_file, "a") as log_f:
            subprocess.run(command, stdout=log_f, stderr=subprocess.STDOUT)
        os.rename(f"{dataset_sep}.col_runs.bv", f"{dataset_sep}.tnl.s10.col_runs.bv")
        os.rename(f"{dataset_sep}.col_ids", f"{dataset_sep}.tnl.s10.col_ids")
        os.rename(f"{dataset_sep}.full.col_ids", f"{dataset_sep}.tnl.s10.full.col_ids")

    dataset_fa_dir = os.getenv("DATA_DIR", "")
    if not os.path.exists(os.path.join(dataset_fa_dir, f"{dataset_sep}.colbwt")):
        log("[bwt_frag]", f"Make Movi output for {filename} with tunnels and split rate 10")
        command = [movi_build_split, "build", "-f", dataset_sep, "-i", os.path.join(dataset_fa_dir, "movi_split_tnl_s10")]
        with open(log_file, "a") as log_f:
            subprocess.run(command, stdout=log_f, stderr=subprocess.STDOUT)

    log("Finished processing", filename)
    

def build(args):
    log_file = os.path.join(os.getenv("DATA_DIR", ""), filename)
    with open(log_file, "w") as f:
        f.write("Processing fasta files\n")

    dataset_sep = f"{filename}.fa"
    dataset_files = "files.txt"

def query(args):
    # Add your query logic here
    print("Running query mode with arguments:", args)

def main():
    parser = argparse.ArgumentParser(description="col-bwt: A co-linear BWT representations supporting chain statistics")
    subparsers = parser.add_subparsers(dest='command', help='Subcommands')

    # Build subcommand
    parser_prepare = subparsers.add_parser('prepare', help='Build the BWT and find multi-MUMs')
    group = parser_prepare.add_mutually_exclusive_group(required=True)
    group.add_argument('input', help='input file name', type=str, nargs='+')
    group.add_argument('-i', '--input-file', help='path to a file-list of genomes to use', type=str)

    parser_prepare.add_argument('-o', '--output', help='output prefix path', type=str, required=True)
    parser_prepare.add_argument('-A', '--arrays-out', help='write the full LCP, BWT, and SA to file', action='store_true')
    parser_prepare.add_argument('-v', help='verbose', action='store_true')
    parser_prepare.set_defaults(func=prepare)

    # Query subcommand
    parser_query = subparsers.add_parser('query', help='Query the BWT')
    parser_query.add_argument('input', help='input file name', type=str)
    parser_query.add_argument('-p', help='pattern to search for', type=str, required=True)
    parser_query.add_argument('-v', help='verbose', action='store_true')
    parser_query.set_defaults(func=query)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
    else:
        args.func(args)

if __name__ == "__main__":
    main()