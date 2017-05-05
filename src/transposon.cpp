// -*- mode: c++; coding: utf-8 -*-
/*! @file transposon.cpp
    @brief Implementation of Transposon class
*/
#include "transposon.hpp"

#include <iostream>

#include "wtl/debug.hpp"

namespace tek {

double Transposon::ALPHA_ = 0.8;
unsigned int Transposon::BETA_ = 24;

std::ostream& Transposon::write(std::ostream& ost) const {
    ost << sequence_;
    return ost;
}

std::ostream& operator<<(std::ostream& ost, const Transposon& ind) {
    return ind.write(ost);
}

void Transposon::unit_test() {HERE;
    Transposon te;
    std::cout << te << std::endl;
}


} // namespace tek
