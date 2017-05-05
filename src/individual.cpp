// -*- mode: c++; coding: utf-8 -*-
/*! @file individual.cpp
    @brief Implementation of Individual class
*/
#include "individual.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

#include "wtl/debug.hpp"

namespace tek {

double Individual::XI_ = 1e-4;

Individual::Individual(const haploid_t& hap0, const haploid_t& hap1)
: chromosomes_{hap0, hap1} {}

double Individual::fitness() const {
    std::valarray<double> s_gpj(NUM_SITES);  // TODO: class variable
    const auto z = genotype();
    const std::valarray<double> v = 1.0 - z * s_gpj;
    double s_gp = 1.0;
    s_gp -= std::accumulate(std::begin(v), std::end(v), 1.0, std::multiplies<double>());
    return (1.0 - s_gp) * (1.0 - s_cn(z.sum()));
}

std::valarray<double> Individual::genotype() const {
    std::valarray<double> z(NUM_SITES);
    for (size_t j=0; j<NUM_SITES; ++j) {
        if (chromosomes_.first[j]) ++z[j];
        if (chromosomes_.second[j]) ++z[j];
    }
    return z;
}

double Individual::s_cn(const unsigned int n) const {
    return std::pow(XI_ * n, TAU_);
}

std::ostream& Individual::write(std::ostream& ost) const {
    ost << fitness();
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
