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
std::uniform_int_distribution<size_t> Haploid::SITES_DIST_(0, Haploid::NUM_SITES - 1U);
std::shared_ptr<Transposon> Haploid::ORIGINAL_TE_ = std::make_shared<Transposon>();

size_t Haploid::random_index() {
    return SITES_DIST_(wtl::sfmt());
}

Haploid::Haploid(): sites_(NUM_SITES) {
    sites_[random_index()] = ORIGINAL_TE_;
}

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
    std::vector<size_t> positions(num_chiasma + 1U);
    for (size_t i=0; i<num_chiasma; ++i) {
        positions[i] = random_index();
    }
    positions[num_chiasma] = -1U;
    std::sort(positions.begin(), positions.end());
    bool flg = false;
    for (size_t i=0, j=0; j<NUM_SITES; ++j) {
        if (j == positions[i]) {
            flg ^= flg;
            ++i;
        }
        if (flg) {
            auto tmp = sites_[j];
            sites_[j] = other.sites_[j];
            other.sites_[j] = tmp;
        }
    }
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
        const bool is_deactivating = bern_indel(wtl::sfmt());
        if (num_mutations > 0 || is_deactivating) {
            p = std::make_shared<Transposon>(*p);
        }
        for (size_t i=0; i<num_mutations; ++i) {
            p->mutate();
        }
        if (is_deactivating) {
            p->indel();
        }
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
