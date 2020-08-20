# Ram
[![Build status for c++/clang++](https://travis-ci.com/dbabojelic/ram.svg?branch=master)](https://travis-ci.org/dbabojelic/ram)


Ram is a c++ implementation of [minimap](https://github.com/lh3/minimap) with few modifications.

# Usage

To build ram run the following commands:
```bash
git clone --recursive https://github.com/dbabojelic/ram.git ram
cd ram && mkdir build && cd build
cmake -Dram_build_executable=ON -DCMAKE_BUILD_TYPE=Release .. && make
./bin/ram
```
which will display the following usage:
```bash
usage: ram [options ...] <target> [<sequences>]

  # default output is stdout
  <target>/<sequences>
    input file in FASTA/FASTQ format (can be compressed with gzip)

  options will be applied sequentially as specified, example:
  $ ram -w10 -k19 -w5 reads.fastq
  will result in w = 5

  options:
    -k, --kmer-length <int>
      default: 15
      length of minimizers
    -w, --window-length <int>
      default: 5
      length of sliding window from which minimizers are found
    -H, --hpc
      Use homopolymer-compressed (HPC) minimizers
    -r, --robust-winnowing
      Use robust winnowing while extracting minimizers (idea taken from Winnowmap)
    -f, --frequency-threshold <float>
      default: 0.001
      threshold for ignoring most frequent minimizers
    -M, --Micromize
      use only a portion of all minimizers
    -p, --Micromize-factor <float>
      Expect to get a floating number between 0 and 1
      When using micromizers reduce the number of minimizers to <float> smallest ones
      If zero: number of taken micromizers will be bounded by the value of sequence_len / k
      default: 0
    -N, --Micromize-extend <int>
      when using micromizers always take first and last <int> minimizers
      default: 0
    -K, --begin-end <int>
      when greater than zero, begin-end strategy will be used
      default: 0
    -m <int>
      default: 100
      discard chains with chaining score less than <int>
    -g <int>
      default: 10000
      stop chain elongation if there are no minimizer withing <int>-BP
    -n <int>
      default: 4
      discard chains consisting of less then <int> minimizers
    -b --best-n <int>
      default: 0
      choose only <int> best hits; if zero all hits will be chosen
    -i, --reduce-win-sz <int>
      default: 0
      if zero does nothing; otherwise one more hierarchical level of minimizing procedure is applied (with given window size)
    -x, --preset-options ava|map
      default: none
      preset options; applies multiple options at the same time;
      this options will be overwritten if used with other options;
      available preset options strings:
          ava: all-vs-all alignment (-k19 -w5 -m100 -g10000 -n4)
          map: read to reference mapping (-k19 -w10 -m40 -g5000 -n3 -b5)
    -t, --threads <int>
      default: 1
      number of threads
    --version
      prints the version number
    -h, --help
      prints the usage
```

If you would like to add ram as a library to your project via CMake, add the following:
```cmake
if (NOT TARGET ram)
  add_subdirectory(<path_to_submodules>/ram EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(<your_exe> ram)
```

#### Dependencies

- gcc 4.8+ or clang 3.5+
- cmake 3.9+
- zlib (for binary only)

## Unit tests

To build ram unit tests run the following commands:

```bash
git clone https://github.com/dbabojelic/ram.git ram
cd ram && mkdir build && cd build
cmake -Dram_build_tests=ON -DCMAKE_BUILD_TYPE=Release .. && make
./bin/ram_test
```

#### Dependencies
- gtest

## Acknowledgement

This work has been supported in part by the European Regional Development Fund under the grant KK.01.1.1.01.0009 (DATACROSS) and in part by the Croatian Science Foundation under the project Single genome and metagenome assembly (IP-2018-01-5886).
