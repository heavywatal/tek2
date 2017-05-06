// -*- mode: c++; coding: utf-8 -*-
/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"

#include <iostream>

#include <wtl/debug.hpp>
#include <wtl/prandom.hpp>

namespace tek {

Population::Population(const size_t n)
: gametes_(2 * n, haploid_t(Individual::NUM_SITES)) {HERE;
}

void Population::step() {HERE;
    std::vector<haploid_t> nextgen;
    nextgen.reserve(gametes_.size());
    std::uniform_int_distribution<size_t> unif(0, gametes_.size() - 1U);
    while (nextgen.size() < gametes_.size()) {
        const auto& egg = gametes_[unif(wtl::sfmt())];
        const auto& sperm = gametes_[unif(wtl::sfmt())];
        Individual diploid(egg, sperm);
        std::bernoulli_distribution bernoulli(diploid.fitness());
        if (!bernoulli(wtl::sfmt())) continue;
        diploid.transpose();
        diploid.recombine();
        diploid.mutate();
        auto genome = diploid.extract_genome();
        nextgen.push_back(std::move(genome.first));
        nextgen.push_back(std::move(genome.second));
    }
    gametes_.swap(nextgen);
}

std::ostream& Population::write(std::ostream& ost) const {HERE;
    for (const auto& x: gametes_) {
        ost << x << std::endl;
    }
    return ost;
}

std::ostream& operator<<(std::ostream& ost, const Population& pop) {
    return pop.write(ost);
}

void Population::unit_test() {HERE;
    Population pop(4);
    std::cout << pop << std::endl;
}

} // namespace tek
