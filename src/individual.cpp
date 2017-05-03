// -*- mode: c++; coding: utf-8 -*-
/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "individual.hpp"

#include <iostream>

#include "wtl/debug.hpp"

namespace tek {

std::ostream& Individual::write(std::ostream& ost) const {
    ost << fitness_;
    return ost;
}

std::ostream& operator<<(std::ostream& ost, const Individual& ind) {
    return ind.write(ost);
}

void Individual::unit_test() {HERE;
    Individual ind;
    std::cout << ind << std::endl;
}


} // namespace tek
