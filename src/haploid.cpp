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
double Haploid::MEAN_SELECTION_COEF_ = 1e-4;
std::valarray<double> Haploid::SELECTION_COEFS_GP_(Haploid::NUM_SITES);
std::uniform_int_distribution<size_t> Haploid::SITES_DIST_(0, Haploid::NUM_SITES - 1U);
std::shared_ptr<Transposon> Haploid::ORIGINAL_TE_ = std::make_shared<Transposon>();

void Haploid::set_SELECTION_COEFS_GP() {
    std::bernoulli_distribution bernoulli(PROP_FUNCTIONAL_SITES_);
    std::exponential_distribution<double> expo_dist(MEAN_SELECTION_COEF_);
    for (size_t i=0; i<NUM_SITES; ++i) {
        if (bernoulli(wtl::sfmt())) {
            SELECTION_COEFS_GP_[i] = expo_dist(wtl::sfmt());
        }
    }
}

size_t Haploid::random_index() {
    return SITES_DIST_(wtl::sfmt());
}

Haploid::Haploid(): sites_(NUM_SITES) {
    sites_[random_index()] = ORIGINAL_TE_;
}

double Haploid::selection_coef_cn(const unsigned int n) {
    return std::pow(XI_ * n, TAU_);
}

double Haploid::fitness(const Haploid& other) const {
    const auto z = valarray() + other.valarray();
    const std::valarray<double> v = 1.0 - z * SELECTION_COEFS_GP_;
    double s_gp = 1.0;
    s_gp -= std::accumulate(std::begin(v), std::end(v), 1.0, std::multiplies<double>());
    return (1.0 - s_gp) * (1.0 - selection_coef_cn(z.sum()));
}

std::vector<std::shared_ptr<Transposon>> Haploid::transpose() {
    std::bernoulli_distribution bern_excision(EXCISION_RATE_);
    std::vector<std::shared_ptr<Transposon>> copying_transposons;
    for (auto& p: sites_) {
        if (!p) continue;
        std::bernoulli_distribution bern(p->transposition_rate());
        if (bern(wtl::sfmt())) {
            copying_transposons.push_back(p);
        }
        if (bern_excision(wtl::sfmt())) {p.reset();}
    }
    return copying_transposons;
}

void Haploid::transpose(Haploid& other) {
    auto copying_transposons = this->transpose();
    {
        auto tmp = other.transpose();
        copying_transposons.insert(copying_transposons.end(),
            std::make_move_iterator(tmp.begin()),
            std::make_move_iterator(tmp.end()));
    }
    std::bernoulli_distribution coin_dist(0.5);
    for (auto& p: copying_transposons) {
        if (coin_dist(wtl::sfmt())) {
            auto& dest = this->sites_[random_index()];
            if (dest) continue;
            dest = std::move(p);
        } else {
            auto& dest = other.sites_[random_index()];
            if (dest) continue;
            dest = std::move(p);
        }
    }
    // TODO: should we choose targets from empty sites?
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
