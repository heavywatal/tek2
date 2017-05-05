// -*- mode: c++; coding: utf-8 -*-
/*! @file transposon.cpp
    @brief Implementation of Transposon class
*/
#include "transposon.hpp"

#include <iostream>

#include "wtl/debug.hpp"
#include "wtl/math.hpp"
#include "wtl/prandom.hpp"

namespace tek {

double Transposon::ALPHA_ = 0.8;
unsigned int Transposon::BETA_ = 24U;

void Transposon::mutate() {
    std::uniform_int_distribution<unsigned int> dist(0U, SEQUENCE_LENGTH_ - 1U);
    sequence_.flip(dist(wtl::sfmt()));
}

double Transposon::activity() const {
    if (has_indel_) return 0.0;
    const double diff = sequence_.count() * OVER_L_;
    const double threshold = 1.0 - ALPHA_;  // TODO: class variable
    if (diff >= threshold) return 0.0;
    return wtl::pow((threshold - diff) / threshold, BETA_);
}

std::ostream& Transposon::write(std::ostream& ost) const {
    return ost << has_indel_ << sequence_;
}

std::ostream& operator<<(std::ostream& ost, const Transposon& ind) {
    return ind.write(ost);
}

void Transposon::unit_test() {HERE;
    Transposon te;
    te.mutate();
    te.mutate();
    std::cout << te << std::endl;
    std::cout << te.activity() << std::endl;
    te.indel();
    std::cout << te << std::endl;
    std::cout << te.activity() << std::endl;
}

} // namespace tek
