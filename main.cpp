// -*- mode: c++; coding: utf-8 -*-
/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include "src/program.hpp"
#include <iostream>
#include <stdexcept>

//! Just instantiate and run Program
int main(int argc, char* argv[]) {
    try {
        tek::Program program(argc, argv);
        program.run();
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}
