// -*- mode: c++; coding: utf-8 -*-
/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"

#include <iostream>

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/prandom.hpp>

namespace tek {

Population::Population(const size_t size, const size_t num_founders) {HERE;
    Haploid founder;
    gametes_.reserve(size * 2U);
    founder.init_founder();
    for (size_t i=0; i<num_founders; ++i) {
        gametes_.push_back(founder);
    }
    gametes_.resize(size * 2U);
}

bool Population::evolve(const size_t max_generations) {HERE;
    for (size_t t=0; t<max_generations; ++t) {
        bool is_recording = (t % 5U > 0U);
        const auto activities = step(is_recording);
        std::cerr << "." << std::flush;
        if (is_extinct()) {
            std::cerr << "Extinction!" << std::endl;
            return false;
        }
    }
    std::cerr << *this << std::endl;
    return true;
}

std::vector<double> Population::step(const bool is_recording) {
    std::vector<Haploid> nextgen;
    nextgen.reserve(gametes_.size());
    std::uniform_int_distribution<size_t> unif(0, gametes_.size() - 1U);
    std::vector<double> activities;
    if (is_recording) activities.reserve(gametes_.size());
    while (nextgen.size() < gametes_.size()) {
        auto egg = gametes_[unif(wtl::sfmt())];
        auto sperm = gametes_[unif(wtl::sfmt())];
        std::bernoulli_distribution bernoulli(egg.fitness(sperm));
        if (!bernoulli(wtl::sfmt())) continue;
        egg.mutate(sperm);
        if (is_recording) {
            egg.push_back_activities(&activities);
            sperm.push_back_activities(&activities);
        }
        nextgen.push_back(egg);
        nextgen.push_back(sperm);
    }
    gametes_.swap(nextgen);
    return activities;
}

bool Population::is_extinct() const {
    for (const auto& x: gametes_) {
        if (x.has_transposon()) return false;
    }
    return true;
}

std::ostream& Population::write(std::ostream& ost) const {HERE;
    return ost << gametes_ << std::endl;
}

std::ostream& operator<<(std::ostream& ost, const Population& pop) {
    return pop.write(ost);
}

void Population::unit_test() {HERE;
    Population pop(6);
    std::cout << pop << std::endl;
    std::cout << pop.step() << std::endl;
    std::cout << pop << std::endl;
}

} // namespace tek
