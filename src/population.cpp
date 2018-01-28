/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"
#include "transposon.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/concurrent.hpp>
#include <wtl/random.hpp>
#include <sfmt.hpp>
#include <json.hpp>

#include <unordered_map>
#include <algorithm>
#include <mutex>

namespace tek {

Population::Population(const size_t size, const size_t num_founders, const unsigned int concurrency)
: concurrency_(concurrency) {HERE;
    Haploid::initialize(size, THETA, RHO);
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
            if (Transposon::can_speciate()) {
                supply_new_species();
            }
            if (static_cast<bool>(flags & Recording::activity)) {
                auto ioflag = (t > record_interval) ? std::ios::app : std::ios::out;
                wtl::ozfstream ozf("activity.tsv.gz", ioflag);
                write_activity(ozf, t, t == record_interval);
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
                std::ostringstream outfile;
                outfile << "generation_" << wtl::setfill0w(5) << t << ".fa.gz";
                wtl::ozfstream ozf(outfile.str());
                for (size_t i=0u; i<20u; ++i) {
                    gametes_[i].write_fasta(ozf);
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
    const size_t num_gametes = gametes_.size();
    static Haploid::URBG seeder(std::random_device{}());
    static wtl::ThreadPool pool(concurrency_);
    static std::vector<Haploid> nextgen;
    nextgen.reserve(num_gametes);
    std::vector<double> fitness_record;
    fitness_record.reserve(num_gametes);
    std::mutex mtx;
    auto task = [&,this](const Haploid::URBG::result_type seed) {
        std::uniform_int_distribution<size_t> dist_idx(0u, num_gametes / 2u - 1u);
        Haploid::URBG engine(seed);
        while (true) {
            const size_t mother_idx = dist_idx(engine);
            size_t father_idx = 0u;
            while ((father_idx = dist_idx(engine)) == mother_idx) {;}
            const auto& mother_lchr = gametes_[2u * mother_idx];
            const auto& mother_rchr = gametes_[2u * mother_idx + 1u];
            const auto& father_lchr = gametes_[2u * father_idx];
            const auto& father_rchr = gametes_[2u * father_idx + 1u];
            auto egg   = mother_lchr.gametogenesis(mother_rchr, engine);
            auto sperm = father_lchr.gametogenesis(father_rchr, engine);
            const double fitness = egg.fitness(sperm);
            if (fitness < wtl::generate_canonical(engine) * previous_max_fitness) continue;
            egg.transpose_mutate(sperm, engine);
            std::lock_guard<std::mutex> lock(mtx);
            if (nextgen.size() >= num_gametes) break;
            fitness_record.push_back(fitness);
            nextgen.push_back(std::move(egg));
            nextgen.push_back(std::move(sperm));
        }
    };
    for (size_t i=0u; i<concurrency_; ++i) {
        pool.submit(std::bind(task, seeder()));
    }
    pool.wait();
    gametes_.swap(nextgen);
    nextgen.clear();
    return fitness_record;
}

void Population::supply_new_species() {
    std::unordered_map<uint_fast32_t, TransposonFamily> counter;
    for (const auto& chr: gametes_) {
        for (const auto& p: chr) {
            counter[p.second->species()] += *p.second;
        }
    }
    std::unordered_map<uint_fast32_t, Transposon> centers;
    for (const auto& p: counter) {
        centers[p.first] = p.second.majority();
    }

    Transposon::INTERACTION_COEFS_clear();
    for (const auto& px: centers) {
        for (const auto& py: centers) {
            if (px.first < py.first) {
                Transposon::INTERACTION_COEFS_emplace(px.first, py.first, px.second * py.second);
            }
        }
    }

    Transposon* farthest = nullptr;
    size_t max_distance = 0;
    for (const auto& chr: gametes_) {
        for (const auto& p: chr) {
            if (p.second->activity() < 0.01) continue;
            size_t distance = (*p.second - centers[p.second->species()]);
            if (distance > max_distance) {
                max_distance = distance;
                farthest = p.second.get();
            }
        }
    }

    if (farthest &&
        std::all_of(centers.begin(), centers.end(), [farthest](const auto& p) {
            return farthest->is_far_enough_from(p.second);
        })) {
        farthest->speciate();
        DCERR(*farthest);
        // NOTE: the center of the original species has not been adjusted
        for (const auto& p: centers) {
            Transposon::INTERACTION_COEFS_emplace(p.first, farthest->species(), p.second * *farthest);
        }
    }
}

bool Population::is_extinct() const {
    return std::all_of(gametes_.begin(), gametes_.end(), [](const Haploid& x) {
        return x.empty();
    });
}

void Population::write_activity(std::ostream& ost, const size_t time, const bool header) const {
    std::map<uint_fast32_t, std::map<double, uint_fast32_t>> counter;
    for (const auto& chr: gametes_) {
        for (const auto& p: chr) {
            const auto& te = *p.second;
            ++counter[te.species()][te.activity()];
        }
    }
    if (header) {
        ost << "generation\tspecies\tactivity\tcopy_number\n";
    }
    for (const auto& sp: counter) {
        auto species = sp.first;
        for (const auto& act_cnt: sp.second) {
            ost << time << "\t" << species << "\t"
                << act_cnt.first << "\t" << act_cnt.second << "\n";
        }
    }
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
        p.first->write_metadata(ost << ">");
        ost << " copy_number=" << p.second << "\n";
        p.first->write_sequence(ost) << "\n";
    }
    return ost;
}

//! shortcut << Population::gametes_
std::ostream& operator<<(std::ostream& ost, const Population& pop) {
    return ost << pop.gametes_;
}

} // namespace tek
