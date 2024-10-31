# col-bwt
```console
               /\                 /\                          /\     
              /**\           /\  /**\                        /**\    
             /****\   /\    /**\/****\                    __/****\     
            /      \ /**\__/   /      \                  /        \ 
           /  /\    /    \    /        \ _       _      /      _   \__
          /  /  \  /      \  /  ___ ___ | |     | |____/     _| |_    \ 
         /  /    \/ /\     \/  / __/ _ \| |_____| '_ \ \ /\ / / __|    \ 
        /  /      \/  \/\   \ | (_| (_) | |_____| |_) \ V  V /| |_      \ 
      _/__/_______/___/__\___\_\___\___/|_|     |_.__/ \_/\_/  \__|______\_

                                                        Alpha Version 0.2.1
```

Pangenomic index computing PMLs and chain staistics by using multi-MUMs. 
# Description
The BWT is constructed to identify *sub-runs* which correspond to the multi-maximal unique matches among the sequences; we then index this BWT using *Movi*[1]. This supports computing *pseudo matching lengths*[2] in O(m) time for a pattern $P[1..m]$. When these sub-runs can be *tunneled*[3] we construct an index in $O(r+n/d)$-space, where $n$ is the length of the text, $r$ is the number of runs of its $BWT$, and $d$ is the number of documents (sequences) in the collection. For pangenome indexes, $\frac{n}{d}$ is the averafe genome length, which does not grow substantially with the number of sequences.

Chain statistics describe which conserved region an exact match falls in; we can chain seeds during query time, linear with respect to the pattern. This results in improved read classification on pangenome indexes when multi-MUM coverage is high. On 64 human haplotypes of HPRC, ~80% of bases belong to a multi-MUM.

Our index uses *col* as it abbreviates co-linearity (not column). A col is also the lowest point on a ridge between two peaks, which describes our index aptly. In construction, we carve runs from a larger range due to their heightened similarity. In application, we bridge gaps between peaks when they appear to be of the same landscape.

# How-to
### Download and Compile
```console
git clone https://github.com/drnatebrown/col-bwt.git
cd col-bwt

mkdir build && cd build
cmake ..
make
```

### Commands
Python3 is required to run the ``col-bwt`` script.
See help for available commands:
```console
./col-bwt -h
```
#### Build Command:
Use ``-h`` for help
```console
./col-bwt build [-h] [-i INPUT] -o OUTPUT [-r] [-m MODE] [-s SUB_SAMPLE] [-l MIN_MUM]
                     [-v] [--force] [--keep] [--clean]
                     [fastas ...]
```
Example:
```console
./col-bwt build -o ./data/index_files -r -m tunnels -s 10 ./data/seq1.fa ./data/seq2.fa
```
Builds a ``col-bwt`` stored in directory ``./data/index_files``, using multi-MUM tunnels and sub-sampling rate 10, over the documents ``seq1.fa``,``seq2.fa`` and their reverse compliment.

Input Format:
``-i`` is an override expecting a path to a file of genomes, one per line, with their space delimited "class", e.g.
```console
/path/to/seq1.fa 1
/path/to/seq2.fa 2
/path/to/seq3.fa 3
/path/to/seq4.fa 4
```
If run with positional arguments, a file ``[PREFIX]_filelist.txt`` is generated. It is highly recommended to use ``-i`` with this file when running the build script again. Since many steps are skipped if intermediate files are already generated, this ensures they are to the expected dataset.

#### Query Command:
Use ``-h`` for help
```console
col-bwt query [-h] -p PATTERN index
```
Example:
```console
./col-bwt query -o -p ./data/pattern.fa ./data/index_files
```
Generates ``pattern.fa.split.pml.bin`` and ``pattern.fa.cid.bin``, the *pseudo matching lengths* and *chain statistics respectively*.

# Thirdparty
\* *denotes fork modified for this tool*
* [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)
    * [divufsort](https://github.com/simongog/libdivsufsort)
* [mumemto](https://github.com/drnatebrown/mumemto.git)\*
* [Movi](https://github.com/drnatebrown/Movi.git)\*

# Citations
[1] Zakeri, M., Brown, N. K., Ahmed, O. Y., Gagie, T., & Langmead, B. (2023). Movi: a fast and cache-efficient full-text pangenome index. bioRxiv.  
[2] Uwe Baier (2018). On Undetected Redundancy in the Burrows-Wheeler Transform. In Annual Symposium on Combinatorial Pattern Matching, CPM 2018, 3:1–3:15.  
[3] Ahmed, O., Rossi, M., Kovaka, S., Schatz, M. C., Gagie, T., Boucher, C., & Langmead, B. (2021). Pan-genomic matching statistics for targeted nanopore sequencing. Iscience, 24(6).  
