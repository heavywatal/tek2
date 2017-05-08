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

void Population::step() {
    std::vector<Haploid> nextgen;
    nextgen.reserve(gametes_.size());
    std::uniform_int_distribution<size_t> unif(0, gametes_.size() - 1U);
    while (nextgen.size() < gametes_.size()) {
        auto egg = gametes_[unif(wtl::sfmt())];
        auto sperm = gametes_[unif(wtl::sfmt())];
        std::bernoulli_distribution bernoulli(egg.fitness(sperm));
        if (!bernoulli(wtl::sfmt())) continue;
        egg.mutate(sperm);
        nextgen.push_back(egg);
        nextgen.push_back(sperm);
    }
    gametes_.swap(nextgen);
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
    pop.step();
    std::cout << pop << std::endl;
}

} // namespace tek
