// -*- mode: c++; coding: utf-8 -*-
/*! @file individual.cpp
    @brief Implementation of Individual class
*/
#include "individual.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

#include "wtl/debug.hpp"
#include "wtl/iostr.hpp"
#include "wtl/prandom.hpp"

namespace tek {

std::ostream& operator<<(std::ostream& ost, const haploid_t& chr) {
    size_t occupied = 0;
    size_t active = 0;
    for (const auto& p: chr) {
        if (p) {
            ++occupied;
            if (p->transposition_rate() > 0.0) ++active;
        }
    }
    return ost << "(" << active << "/" << occupied << ")";
}

double Individual::XI_ = 1e-4;
double Individual::EXCISION_RATE_ = 1e-6;

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
        if (genome_.first[j]) ++z[j];
        if (genome_.second[j]) ++z[j];
    }
    return z;
}

double Individual::s_cn(const unsigned int n) const {
    return std::pow(XI_ * n, TAU_);
}

void Individual::transpose() {
    std::bernoulli_distribution bern_excision(EXCISION_RATE_);
    for (auto& p: genome_.first) {
        if (!p) continue;
        std::bernoulli_distribution bern(p->transposition_rate());
        if (bern(wtl::sfmt())) {
            // TODO: record the transposon and its new site
        }
        if (bern_excision(wtl::sfmt())) {p.reset();}
    }
    // TODO: second
    // TODO: actual transposition
}

void Individual::recombine() {
    // TODO: population parameters
    constexpr double rho = 200;
    constexpr size_t popsize = 1000;
    const double c = rho / popsize / 4.0;
    std::poisson_distribution<size_t> poisson(c * NUM_SITES); // -1?
    const size_t num_chiasma = poisson(wtl::sfmt());
    if (num_chiasma == 0U) return;
    std::uniform_int_distribution<size_t> unif(0, NUM_SITES - 1);
    std::vector<size_t> positions(num_chiasma);
    for (size_t i=0; i<num_chiasma; ++i) {
        positions[i] = unif(wtl::sfmt());
    }
    std::sort(positions.begin(), positions.end());
    // TODO: make new haploid_t instances
}

void Individual::mutate() {
    // TODO: population parameters
    constexpr double theta = 0.01;
    constexpr size_t popsize = 1000;
    const double mu = theta / popsize / 4.0;  // per TE site?
    std::poisson_distribution<size_t> poisson(mu);
    std::bernoulli_distribution bern_indel(mu * INDEL_RATIO_);
    for (auto& p: genome_.first) {
        if (!p) continue;
        const size_t num_mutations = poisson(wtl::sfmt());
        for (size_t i=0; i<num_mutations; ++i) {
            p->mutate();
        }
        if (bern_indel(wtl::sfmt())) {p->indel();}
    }
    // TODO: second
}

std::ostream& Individual::write(std::ostream& ost) const {
    ost << genome_;
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
