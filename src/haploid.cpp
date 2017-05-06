// -*- mode: c++; coding: utf-8 -*-
/*! @file haploid.cpp
    @brief Implementation of Haploid class
*/
#include "haploid.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

#include "wtl/debug.hpp"
#include "wtl/iostr.hpp"
#include "wtl/prandom.hpp"

namespace tek {

double Haploid::XI_ = 1e-4;
double Haploid::EXCISION_RATE_ = 1e-6;

Haploid::Haploid(): sites_(NUM_SITES) {}

double Haploid::s_cn(const unsigned int n) {
    return std::pow(XI_ * n, TAU_);
}

double Haploid::fitness(const Haploid& other) const {
    std::valarray<double> s_gpj(NUM_SITES);  // TODO: class variable
    const auto z = valarray() + other.valarray();
    const std::valarray<double> v = 1.0 - z * s_gpj;
    double s_gp = 1.0;
    s_gp -= std::accumulate(std::begin(v), std::end(v), 1.0, std::multiplies<double>());
    return (1.0 - s_gp) * (1.0 - s_cn(z.sum()));
}

void Haploid::transpose(Haploid& other) {
    std::bernoulli_distribution bern_excision(EXCISION_RATE_);
    for (auto& p: sites_) {
        if (!p) continue;
        std::bernoulli_distribution bern(p->transposition_rate());
        if (bern(wtl::sfmt())) {
            // TODO: record the transposon and its new site
        }
        if (bern_excision(wtl::sfmt())) {p.reset();}
    }
    // TODO: other
    // TODO: actual transposition
}

void Haploid::recombine(Haploid& other) {
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

void Haploid::mutate() {
    // TODO: population parameters
    constexpr double theta = 0.01;
    constexpr size_t popsize = 1000;
    const double mu = theta / popsize / 4.0;  // per TE site?
    std::poisson_distribution<size_t> poisson(mu);
    std::bernoulli_distribution bern_indel(mu * INDEL_RATIO_);
    for (auto& p: sites_) {
        if (!p) continue;
        const size_t num_mutations = poisson(wtl::sfmt());
        for (size_t i=0; i<num_mutations; ++i) {
            p->mutate();
        }
        if (bern_indel(wtl::sfmt())) {p->indel();}
    }
}

std::ostream& Haploid::write(std::ostream& ost) const {
    size_t occupied = 0;
    size_t active = 0;
    for (const auto& p: sites_) {
        if (p) {
            ++occupied;
            if (p->transposition_rate() > 0.0) ++active;
        }
    }
    return ost << "(" << active << "/" << occupied << ")";
}

std::ostream& operator<<(std::ostream& ost, const Haploid& x) {
    return x.write(ost);
}

void Haploid::unit_test() {HERE;
    Haploid x;
    std::cout << x << std::endl;
}

} // namespace tek
