#!/usr/bin/env zsh

PATTERN=/insert/path/to/pattern/fasta
INPUT_DIR=insert/input/path
INPUT_NAME=input

MOVI_DIR=/path/to/Movi

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

movi_build_split=$MOVI_DIR/build/movi-split

log_file=$PATTERN.log
echo "Processing pattern files" > $log_file

$dataset_sep = $INPUT_DIR/$INPUT_NAME.fa
      
big_log "[[ Processing $PATTERN ]]"

if [[ ! -d $dataset_sep.colbwt ]]; then
    log "[error]" "Need to run build.sh"
    exit 1
fi
log "[col-bwt]" "Generating PMLs and CIDs using index $INPUT_NAME"
exec $movi_build_split query -i $dataset_sep.colbwt -r $PATTERN >> $log_file

log "Finished processing pattern $PATTERN"