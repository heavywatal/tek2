/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include "program.hpp"
#include <iostream>

//! Just instantiate and run Program
int main(int argc, char* argv[]) {
    try {
        tek::Program program(argc, argv);
        program.run();
    } catch (const tek::exit_success& e) {
        return 0;
    } catch (const std::runtime_error& e) {
        std::cerr << "\nruntime_error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
