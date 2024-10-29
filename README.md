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

# Requirements
* ``python3``

# Thirdparty
\* *denotes fork modified for this tool*
* [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)
    * [divufsort](https://github.com/simongog/libdivsufsort)
* [klib](https://github.com/attractivechaos/klib.git)
* [mumemto](https://github.com/drnatebrown/mumemto.git)\*
* [Movi](https://github.com/drnatebrown/Movi.git)\*
