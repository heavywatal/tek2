// -*- mode: c++; coding: utf-8 -*-
/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"

#include <iostream>
#include <algorithm>

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/prandom.hpp>
#include <wtl/zfstream.hpp>

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
    const size_t record_interval = gametes_.size() * 10U;
    auto oss = wtl::make_oss();
    oss << "time\tactivity\tcopies\n";
    for (size_t t=1; t<=max_generations; ++t) {
        bool is_recording = ((t % record_interval) == 0U);
        const auto counter = step(is_recording);
        if (is_recording) {
            for (const auto& p: counter) {
                oss << t << "\t" << p.first << "\t" << p.second << "\n";
            }
        }
        std::cerr << "." << std::flush;
        if (is_extinct()) {
            std::cerr << "Extinction!" << std::endl;
            return false;
        }
    }
    wtl::ozfstream("activities.tsv") << oss.str();
    std::cerr << *this << std::endl;
    return true;
}

std::map<double, unsigned int> Population::step(const bool is_recording) {
    std::vector<Haploid> nextgen;
    nextgen.reserve(gametes_.size());
    std::uniform_int_distribution<size_t> unif(0, gametes_.size() - 1U);
    std::map<double, unsigned int> counter;
    while (nextgen.size() < gametes_.size()) {
        auto egg = gametes_[unif(wtl::sfmt())];
        auto sperm = gametes_[unif(wtl::sfmt())];
        std::bernoulli_distribution bernoulli(egg.fitness(sperm));
        if (!bernoulli(wtl::sfmt())) continue;
        egg.mutate(sperm);
        if (is_recording) {
            egg.count_activities(&counter);
            sperm.count_activities(&counter);
        }
        nextgen.push_back(egg);
        nextgen.push_back(sperm);
    }
    gametes_.swap(nextgen);
    return counter;
}

bool Population::is_extinct() const {
    for (const auto& x: gametes_) {
        if (x.has_transposon()) return false;
    }
    return true;
}

void Population::sample() const {
    const size_t sample_size = std::max(gametes_.size() / 100UL, 2UL);
    std::ostringstream oss;
    for (size_t i=0; i<sample_size; ++i) {
        gametes_[i].write_sample(oss);
    }
    std::cerr << oss.str();
}

std::ostream& Population::write(std::ostream& ost) const {HERE;
    return ost << gametes_ << std::endl;
}

std::ostream& operator<<(std::ostream& ost, const Population& pop) {
    return pop.write(ost);
}

void Population::unit_test() {HERE;
    Population pop(6, 6);
    std::cout << pop;
    std::cout << pop.step() << std::endl;
    std::cout << pop;
    pop.sample();
}

} // namespace tek
