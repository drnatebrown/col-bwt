#!/usr/bin/env zsh

# USER DEFINED VARIABLES
DATA_DIR=/home/nbrown99/vast/col_stats/data/hprc
DATA_NAMES=(4 8 16 32 64 88)
CHR_DIR=/home/nbrown99/data/public/vshiv/hprc_scaffolds
CHR_NAMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)

MUMEMTO_DIR=/home/nbrown99/repos/mumemto
MIN_MUM=20

PFP_THRESHOLDS_DIR=/home/nbrown99/repos/pfp-thresholds

BWT_FRAG_DIR=/home/nbrown99/repos/bwt_frag

# DEFINITIONS
# alias Time='/usr/bin/time --format="Wall Time: %e\nMax Memory: %M"'
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

col_directory=$DATA_DIR/col_data
mumemto_prg=$MUMEMTO_DIR/build/mumemto
pfp_thr_prg=$PFP_THRESHOLDS_DIR/build/test/src/pfp_thresholds
col_split_build=$BWT_FRAG_DIR/temp/src/build_FL
col_split_prg=$BWT_FRAG_DIR/temp/src/col_split

# Loop over each file
for chr in $CHR_NAMES; do
    log_file=$DATA_DIR/$chr.tnl.log
    echo "Processing $chr" > $log_file
    for num in $DATA_NAMES; do
        filename=hg.$num.$chr
        dataset=$DATA_DIR/$filename.fa
        dataset_fa_dir=$col_directory/$filename
        dataset_sep=$dataset_fa_dir/$filename.fna

        big_log "[[ Processing $dataset ]]"

        if [[ ! -d $dataset_fa_dir ]]; then
            log "[bash]" "Creating directory $dataset_fa_dir"
            mkdir -p "$dataset_fa_dir"
            # Create the text file for each
            ls -1 $CHR_DIR/$chr/*.fa | head -n $num | awk -v path="$CHR_DIR/$chr" '{print $0 " " NR}' > $dataset_fa_dir/files.txt
        fi

        N=$(wc -l $dataset_fa_dir/files.txt | awk '{print $1}')

        log "Number of genomes: $N"

        mum_file=$dataset_sep.mums
        # if [[ ! -e $mum_file ]]; then
            mumemto_filename=$dataset_fa_dir/$filename
            log "[mumemto]" "Building RLBWT and finding MUMs of length atleast $MIN_MUM on $dataset"
            [ -e "$dataset_sep.parse" ] && pfp_flag="-s" || pfp_flag=""
            [ ! -e "$dataset_sep.bwt.heads" ] && bwt_flag="-R" || bwt_flag=""
            [ ! -e "$dataset_sep.thr" ] && thr_flag="-T" || thr_flag=""
            exec $mumemto_prg mum -i $dataset_fa_dir/files.txt -o $mumemto_filename -K -l $MIN_MUM $pfp_flag $bwt_flag $thr_flag
            mv $mumemto_filename.mums $mum_file
        # fi

        if [[ ! -e $dataset_sep.bwt.heads || ! -e $dataset_sep.thr ]]; then
            log "[pfp-thresholds]" "Writing RLBWT, SA, and thresholds of $dataset"
            exec $pfp_thr_prg -r -f $dataset_fa_dir/$filename.fna
        fi

        if [[ ! -e $dataset_sep.FL_table ]]; then
            log "[bwt_frag]" "Building FL table for $dataset"
            exec $col_split_build $dataset_sep
        fi

        log "[bwt_frag]" "Splitting RLBWT to create runs for MUMs of $dataset"
        exec $col_split_prg $dataset_sep -N $N >> $log_file
    done
done