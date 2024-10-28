#!/usr/bin/env zsh
DATA_DIR=/path/to/data
# Names of fasta in data directory
DATA_NAMES=(seq1.fa seq2.fa) 
OUTPUT_DIR=/path/to/out
OUTPUT_NAME=outname

MUMEMTO_DIR=/path/to/mumemto
MOVI_DIR=/path/to/Movi
COL_BWT_DIR=../

MIN_MUM=20

Time() {
    /usr/bin/time --format="Wall Time: %e\nMax Memory: %M" "$@"
}

exec() {
    echo -e "\e[1;33m>> $@\e[0m"
    Time "$@"
    echo -e "\e[1;33m<<\e[0m"
}

big_log() {
    echo -e "\e[1;31m$1\e[0m"
}

log() {
    if [ $# -eq 2 ]; then
        echo -e "\e[1;32m$1 \e[1;35m$2\e[0m"
    else
        echo -e "\e[1;35m$1\e[0m"
    fi
}

mumemto_prg=$MUMEMTO_DIR/build/mumemto
col_split_build=$COL_BWT_DIR/dev_build/src/build_FL
col_split_prg=$COL_BWT_DIR/dev_build/src/col_split
rlbwt_to_bwt_prg=$COL_BWT_DIR/dev_build/src/to_run_output
movi_build_split=$MOVI_DIR/build/movi-split

log_file=$OUTPUT_DIR/$OUTPUT_NAME.log
echo "Processing fasta files" > $log_file

$filename = $DATA_DIR/$OUTPUT_NAME

if [[ ! -d $OUTPUT_DIR ]]; then
    log "[bash]" "Creating directory $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

dataset_files=$OUTPUT_DIR/files.txt
if [[ ! -e $dataset_files ]]; then
    log "[bash]" "Creating files.txt for $OUTPUT_NAME"
    i=1
    for seq in $DATA_NAMES; do
        echo "$DATA_DIR/$seq $i" >> $dataset_files
        i=$((i + 1))
    done
fi
      
big_log "[[ Processing $OUTPUT_NAME ]]"

N=$(wc -l $dataset_files | awk '{print $1}')

log "Number of sequences: $N"

dataset=$OUTPUT_DIR/$OUTPUT_NAME
dataset_sep=$OUTPUT_DIR/$OUTPUT_NAME.fa
if [[ ! -e $dataset_sep.col_mums ]]; then
    log "[mumemto]" "Building RLBWT, Thresholds; finding MUMs of length atleast $MIN_MUM on $filename"
    [ -e "$dataset_sep.parse" ] && pfp_flag="-s" || pfp_flag=""
    [ ! -e "$dataset_sep.bwt.heads" ] && bwt_flag="-R" || bwt_flag=""
    [ ! -e "$dataset_sep.thr" ] && thr_flag="-T" || thr_flag=""
    exec $mumemto_prg mum -i $dataset_files -o $dataset -r -K -l $MIN_MUM $pfp_flag $bwt_flag $thr_flag >> $log_file
fi

if [[ ! -e $dataset_sep.FL_table ]]; then
    log "[col-bwt]" "Building FL table for $OUTPUT_NAME"
    exec $col_split_build $dataset_sep >> $log_file
fi

if [[ ! -e $dataset_sep.bwt ]]; then
    log "[col-bwt]" "Building raw BWT for $OUTPUT_NAME"
    exec $rlbwt_to_bwt_prg $dataset_sep >> $log_file
fi

if [[ ! -e $dataset_sep.col_runs ]]; then
    log "[col-bwt]" "Splitting RLBWT to create runs for MUMs of $filename with tunnels and split rate 10"
    exec $col_split_prg $dataset_sep -N $N -m tunnels -v -s 10 >> $log_file
fi
if [[ ! -d $dataset_fa_dir/$dataset_sep.colbwt ]]; then
    log "[col-bwt]" "Make Movi output for $filename with tunnels and split rate 10"
    exec $movi_build_split build -f $dataset_sep -i $dataset_fa_dir/movi_split_tnl_s10
fi

log "Finished processing $filename"
