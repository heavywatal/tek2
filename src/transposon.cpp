// -*- mode: c++; coding: utf-8 -*-
/*! @file transposon.cpp
    @brief Implementation of Transposon class
*/
#include "transposon.hpp"

#include <iostream>

#include "wtl/debug.hpp"
#include "wtl/math.hpp"

namespace tek {

double Transposon::ALPHA_ = 0.8;
unsigned int Transposon::BETA_ = 24;

void Transposon::mutate() {
    sequence_.flip(0);
}

double Transposon::activity() const {
    const double d = sequence_.count() * OVER_L_;
    const double a = 1.0 - ALPHA_;  // TODO: class variable
    return wtl::pow((a - d) / a, BETA_);
}

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
    std::cout << te.activity() << std::endl;
    te.mutate();
    std::cout << te << std::endl;
    std::cout << te.activity() << std::endl;
}


} // namespace tek
