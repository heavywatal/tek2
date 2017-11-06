/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"
#include "haploid.hpp"
#include "transposon.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/filesystem.hpp>
#include <wtl/concurrent.hpp>
#include <wtl/prandom.hpp>
#include <sfmt.hpp>

#include <json.hpp>

#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <mutex>

namespace tek {

Population::Population(const size_t size, const size_t num_founders, const unsigned int concurrency)
: concurrency_(concurrency) {HERE;
    gametes_.reserve(size * 2u);
    for (size_t i=0u; i<num_founders; ++i) {
        gametes_.push_back(Haploid::copy_founder());
    }
    gametes_.resize(size * 2u);
}

bool Population::evolve(const size_t max_generations, const size_t record_interval, const Recording flags) {HERE;
    constexpr double margin = 0.1;
    double max_fitness = 1.0;
    for (size_t t=1; t<=max_generations; ++t) {
        bool is_recording = ((t % record_interval) == 0u);
        const auto fitness_record = step(max_fitness);
        max_fitness = *std::max_element(fitness_record.begin(), fitness_record.end());
        max_fitness = std::min(max_fitness + margin, 1.0);
        if (is_recording) {
            std::cerr << "*" << std::flush;
            if (static_cast<bool>(flags & Recording::activity)) {
                auto ioflag = (t > record_interval) ? std::ios::app : std::ios::out;
                wtl::ozfstream ozf("activity.tsv.gz", ioflag);
                if (t == record_interval) {
                    ozf << "generation\tgamete\tactivity\tcopy_number\n";
                }
                for (size_t i=0u; i<gametes_.size(); ++i) {
                    for (const auto& p: gametes_[i].count_activity()) {
                        ozf << t << "\t" << i << "\t"
                                      << p.first << "\t" << p.second << "\n";
                    }
                }
            }
            if (static_cast<bool>(flags & Recording::fitness)) {
                auto ioflag = (t > record_interval) ? std::ios::app : std::ios::out;
                wtl::ozfstream ozf("fitness.tsv.gz", ioflag);
                if (t == record_interval) {
                    ozf << "generation\tfitness\n";
                }
                for (const double w: fitness_record) {
                    ozf << t << "\t" << w << "\n";
                }
            }
            if (static_cast<bool>(flags & Recording::sequence)) {
                std::ostringstream outdir;
                outdir << "generation_" << t;
                wtl::ChDir cd(outdir.str(), true);
                for (size_t i=0u; i<10u; ++i) {
                    std::ostringstream outfile;
                    outfile << "individual_" << i << ".fa.gz";
                    wtl::ozfstream ozf(outfile.str());
                    write_individual(ozf, i);
                }
            }
        } else {
            DCERR("." << std::flush);
        }
        if (is_extinct()) {
            std::cerr << "Extinction!" << std::endl;
            return false;
        }
    }
    std::cerr << std::endl;
    return true;
}

std::vector<double> Population::step(const double previous_max_fitness) {
    static wtl::ThreadPool pool(concurrency_);
    const size_t num_gametes = gametes_.size();
    std::vector<size_t> indices(num_gametes / 2u);
    std::iota(indices.begin(), indices.end(), 0u);
    std::vector<Haploid> nextgen;
    nextgen.reserve(num_gametes);
    std::vector<double> fitness_record;
    fitness_record.reserve(num_gametes);
    std::mutex mtx;
    auto task = [&,this](const size_t seed) {
        wtl::sfmt19937 rng(seed);
        while (true) {
            const auto parents = wtl::sample(indices, 2u, rng);
            const auto& mother_lchr = gametes_[2u * parents[0u]];
            const auto& mother_rchr = gametes_[2u * parents[0u] + 1u];
            const auto& father_lchr = gametes_[2u * parents[1u]];
            const auto& father_rchr = gametes_[2u * parents[1u] + 1u];
            auto egg   = mother_lchr.gametogenesis(mother_rchr, rng);
            auto sperm = father_lchr.gametogenesis(father_rchr, rng);
            const double fitness = egg.fitness(sperm);
            if (fitness < rng.canonical() * previous_max_fitness) continue;
            egg.transpose_mutate(sperm, rng);
            std::lock_guard<std::mutex> lock(mtx);
            if (nextgen.size() >= num_gametes) break;
            fitness_record.push_back(fitness);
            nextgen.push_back(std::move(egg));
            nextgen.push_back(std::move(sperm));
        }
    };
    for (size_t i=0u; i<concurrency_; ++i) {
        pool.submit(std::bind(task, wtl::sfmt()()));
    }
    pool.wait();
    gametes_.swap(nextgen);
    return fitness_record;
}

bool Population::is_extinct() const {
    for (const auto& x: gametes_) {
        if (x.has_transposon()) return false;
    }
    return true;
}

std::ostream& Population::write_summary(std::ostream& ost) const {HERE;
    nlohmann::json record;
    for (const auto& x: gametes_) {
        record.push_back(x.summarize());
    }
    return ost << record;
}

std::ostream& Population::write_fasta(std::ostream& ost) const {HERE;
    std::unordered_map<Transposon*, unsigned int> counter;
    for (const auto& chr: gametes_) {
        for (const auto& p: chr) {
            ++counter[p.second.get()];
        }
    }
    for (const auto& p: counter) {
        p.first->write_fasta(ost, p.second);
    }
    return ost;
}

std::ostream& Population::write_individual(std::ostream& ost, const size_t i) const {
    for (const size_t x: {2u * i, 2u * i + 1u}) {
        for (const auto& p: gametes_.at(x)) {
            p.second->write_fasta(ost);
        }
    }
    return ost;
}

//! shortcut << Population::gametes_
std::ostream& operator<<(std::ostream& ost, const Population& pop) {
    return ost << pop.gametes_;
}

void Population::test() {HERE;
    Population pop(6, 6);
    std::cout << pop.gametes_ << std::endl;
    pop.step();
    std::cout << pop.gametes_ << std::endl;
    pop.write_individual(std::cout, 2);
}

} // namespace tek
