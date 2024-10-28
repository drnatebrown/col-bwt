# col-BWT
Pre-release Version 0.1.0
Full Build/Query Suite WIP

# How-to
### Download and Compile

Compile mumemto:
```console
git clone https://github.com/drnatebrown/mumemto.git
cd mumemto

mkdir build 
cd build && cmake ..
make install
```
Compile Movi:
```console
git clone https://github.com/drnatebrown/Movi.git
cd Movi

mkdir build
cd build
cmake ..
make
```

Compile col-bwt:
```console
git clone https://github.com/drnatebrown/r-index-f.git
cd r-index-f

mkdir build && cd build
cmake ..
make
```

### Build
Update paths/config in ./scripts/build.sh
```console
DATA_DIR=/path/to/data
# Names of fasta in data directory
DATA_NAMES=(seq1.fa seq2.fa) 
OUTPUT_DIR=/path/to/out
OUTPUT_NAME=outname

MUMEMTO_DIR=/path/to/mumemto
MOVI_DIR=/path/to/Movi
```

### Query
Update paths/config in ./scripts/query.sh
```console
PATTERN=/insert/path/to/pattern/fasta
INPUT_DIR=insert/input/path
INPUT_NAME=input

MOVI_DIR=/path/to/Movi
```
