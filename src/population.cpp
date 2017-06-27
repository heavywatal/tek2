// -*- mode: c++; coding: utf-8 -*-
/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"
#include "haploid.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/prandom.hpp>
#include <sfmt.hpp>

#include <iostream>
#include <algorithm>
#include <mutex>
#include <future>

namespace tek {

Population::Population(const size_t size, const size_t num_founders, const unsigned int concurrency)
: concurrency_(concurrency) {HERE;
    gametes_.reserve(size * 2U);
    for (size_t i=0; i<num_founders; ++i) {
        gametes_.push_back(Haploid::copy_founder());
    }
    gametes_.resize(size * 2U);
}

bool Population::evolve(const size_t max_generations) {HERE;
    const size_t record_interval = gametes_.size() / 20U;
    auto oss = wtl::make_oss();
    oss << "time\tactivity\tcopies\n";
    for (size_t t=1; t<=max_generations; ++t) {
        bool is_recording = ((t % record_interval) == 0U);
        const auto counter = step(is_recording);
        if (is_recording) {
            std::cerr << "*" << std::flush;
            for (const auto& p: counter) {
                oss << t << "\t" << p.first << "\t" << p.second << "\n";
            }
        } else {
            std::cerr << "." << std::flush;
        }
        if (is_extinct()) {
            std::cerr << "Extinction!" << std::endl;
            return false;
        }
    }
    wtl::ozfstream("activities.tsv") << oss.str();
    wtl::ozfstream("individuals.tsv.gz") << *this;
    return true;
}

std::map<double, unsigned int> Population::step(const bool is_recording) {
    const size_t num_gametes = gametes_.size();
    std::vector<size_t> indices(num_gametes / 2U);
    std::iota(indices.begin(), indices.end(), 0U);
    std::vector<Haploid> nextgen;
    nextgen.reserve(num_gametes);
    std::map<double, unsigned int> counter;
    std::mutex mtx;
    auto task = [&,this](const size_t seed) {
        wtl::sfmt19937 rng(seed);
        while (true) {
            const auto parents = wtl::sample(indices, 2U, rng);
            const auto& mother_lchr = gametes_[2U * parents[0U]];
            const auto& mother_rchr = gametes_[2U * parents[0U] + 1U];
            const auto& father_lchr = gametes_[2U * parents[1U]];
            const auto& father_rchr = gametes_[2U * parents[1U] + 1U];
            auto egg   = mother_lchr.gametogenesis(mother_rchr, rng);
            auto sperm = father_lchr.gametogenesis(father_rchr, rng);
            const double fitness = egg.fitness(sperm);
            if (fitness < rng.canonical()) continue;
            if (is_recording) {
                std::lock_guard<std::mutex> lock(mtx);
                egg.count_activities(&counter);
                sperm.count_activities(&counter);
            }
            egg.transpose_mutate(sperm, rng);
            std::lock_guard<std::mutex> lock(mtx);
            if (nextgen.size() >= num_gametes) break;
            nextgen.push_back(std::move(egg));
            nextgen.push_back(std::move(sperm));
        }
    };
    std::vector<std::future<void>> futures;
    futures.reserve(concurrency_);
    for (size_t i=0; i<concurrency_; ++i) {
        futures.push_back(std::async(std::launch::async, task, wtl::sfmt()()));
    }
    for (auto& f: futures) f.wait();
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

void Population::test() {HERE;
    Population pop(6, 6);
    std::cout << pop;
    std::cout << pop.step() << std::endl;
    std::cout << pop;
    pop.sample();
}

} // namespace tek
