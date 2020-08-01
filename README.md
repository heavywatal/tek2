# tek2

[![Build status](https://github.com/heavywatal/tek2/workflows/build/badge.svg)](https://github.com/heavywatal/tek2/actions)

*tek2* is a population genetic simulation framework of transposons.


## Reference

W. M. Iwasaki, T. E. Kijima, and H. Innan
(2020) *Mol. Biol. Evol.* **37**(2):355–364.
[pmid:31580443](https://www.ncbi.nlm.nih.gov/pubmed/31580443).
[doi:10.1534/genetics.113.150292](https://doi.org/10.1093/molbev/msz220).
Population Genetics and Molecular Evolution of DNA Sequences in Transposable Elements. II. Accumulation of Variation and Evolution of a New Subfamily.

T. E. Kijima and H. Innan
(2013) *Genetics* **195**(3):957–967.
[pmid:24002643](https://www.ncbi.nlm.nih.gov/pubmed/24002643).
[doi:10.1534/genetics.113.150292](http://dx.doi.org/10.1534/genetics.113.150292).
Population genetics and molecular evolution of DNA sequences in transposable elements. I. A simulation framework.

[Project page on GitHub](https://github.com/heavywatal/tek2)


## Requirements

- Unix-like environment (macOS, Linux, WSL, MinGW on MSYS2, etc.)
- C++14 compiler (clang++ >= Apple LLVM 8.1, g++ >= 5.3)
- [CMake](https://cmake.org/) (>= 3.8.0)

The following libraries are optional or automatically installed:

- [clippson](https://github.com/heavywatal/clippson)
- [cxxwtl](https://github.com/heavywatal/cxxwtl)
- [sfmt-class](https://github.com/heavywatal/sfmt-class)
- [zlib](https://zlib.net)


## Installation

The easiest way is to use [Homebrew](https://brew.sh/).
The following command installs tek2 and all the dependencies:
```sh
brew install heavywatal/tap/tek2
```

Alternatively, you can get the source code from GitHub manually:
```sh
git clone https://github.com/heavywatal/tek2.git
cd tek2/
mkdir build
cd build/
YOUR_PREFIX=${HOME}/local  # or /usr/local
cmake -DCMAKE_INSTALL_PREFIX=$YOUR_PREFIX ..
make -j2
make install
```


## API Document

- [Online documentation generated with doxygen](https://heavywatal.github.io/tek2/)
- @ref params
