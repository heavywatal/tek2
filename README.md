# tek2

[![Build Status](https://travis-ci.com/heavywatal/tek2.svg?branch=master)](https://travis-ci.com/heavywatal/tek2)

T. E. Kijima and H. Innan
(2013) *Genetics* **195**, 3 pp.957-67
[[pmid:24002643](http://www.ncbi.nlm.nih.gov/pubmed/24002643)]
[[doi:10.1534/genetics.113.150292](http://dx.doi.org/10.1534/genetics.113.150292)]
Population genetics and molecular evolution of DNA sequences in transposable elements. I. A simulation framework

[Project page on GitHub](https://github.com/heavywatal/tek2)

@ref params


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

The easiest way is to use [Homebrew](https://brew.sh/)/[Linuxbrew](http://linuxbrew.sh/).
The following command installs tekka and all the dependencies:
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
