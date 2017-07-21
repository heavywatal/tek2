// -*- mode: c++; coding: utf-8 -*-
/*! @file haploid.cpp
    @brief Implementation of Haploid class
*/
#include "haploid.hpp"
#include "transposon.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/prandom.hpp>
#include <sfmt.hpp>

#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>

namespace tek {

double Haploid::XI_ = 1e-4;
double Haploid::EXCISION_RATE_ = 1e-6;
double Haploid::MEAN_SELECTION_COEF_ = 1e-4;

double Haploid::RECOMBINATION_RATE_ = 0.0;
double Haploid::INDEL_RATE_ = 0.0;
std::valarray<double> Haploid::SELECTION_COEFS_GP_(Haploid::NUM_SITES);
std::uniform_int_distribution<uint_fast32_t> Haploid::UNIFORM_SITES_(0, Haploid::NUM_SITES - 1U);
std::poisson_distribution<uint_fast32_t> Haploid::NUM_MUTATIONS_DIST_(0.0);
std::shared_ptr<Transposon> Haploid::ORIGINAL_TE_ = std::make_shared<Transposon>();

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// static functions

namespace po = boost::program_options;

po::options_description Haploid::options_desc() {HERE;
    po::options_description description("Haploid");
    description.add_options()
      ("xi", po::value(&XI_)->default_value(XI_))
      ("nu", po::value(&EXCISION_RATE_)->default_value(EXCISION_RATE_))
      ("lambda", po::value(&MEAN_SELECTION_COEF_)->default_value(MEAN_SELECTION_COEF_));
    return description;
}

void Haploid::set_SELECTION_COEFS_GP() {HERE;
    std::vector<uint_fast32_t> sites(NUM_SITES);
    std::iota(sites.begin(), sites.end(), 0);
    const uint_fast32_t n = NUM_SITES * PROP_FUNCTIONAL_SITES_;
    const auto functional_sites = wtl::sample(sites, n, wtl::sfmt());
    std::exponential_distribution<double> expo_dist(1.0 / MEAN_SELECTION_COEF_);
    for (const auto i: functional_sites) {
        SELECTION_COEFS_GP_[i] = expo_dist(wtl::sfmt());
    }
}

void Haploid::set_parameters(const size_t popsize, const double theta, const double rho) {HERE;
    const double four_n = 4.0 * popsize;
    RECOMBINATION_RATE_ = rho / four_n;
    const double mu = Transposon::LENGTH * theta / four_n;
    INDEL_RATE_ = mu * INDEL_RATIO_;
    NUM_MUTATIONS_DIST_.param(decltype(NUM_MUTATIONS_DIST_)::param_type(mu));
    set_SELECTION_COEFS_GP();
}

Haploid Haploid::copy_founder() {
    // TODO: avoid functional site?
    static auto idx = UNIFORM_SITES_(wtl::sfmt());
    Haploid founder;
    founder.sites_[idx] = ORIGINAL_TE_;
    founder.evaluate_sites();
    return founder;
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

Haploid Haploid::gametogenesis(const Haploid& other, URNG& rng) const {
    Haploid lhalf(*this), rhalf(other);
    lhalf.recombine(rhalf, rng);
    if (rng.canonical() < 0.5) {
        lhalf.evaluate_sites();
        return lhalf;
    } else {
        rhalf.evaluate_sites();
        return rhalf;
    }
}

std::vector<std::shared_ptr<Transposon>> Haploid::transpose(URNG& rng) {
    std::vector<std::shared_ptr<Transposon>> copying_transposons;
    for (auto& p: sites_) {
        if (rng.canonical() < p.second->transposition_rate()) {
            copying_transposons.push_back(p.second);
        }
        if (rng.canonical() < EXCISION_RATE_) {
            sites_.erase(p.first);
        }
    }
    return copying_transposons;
}

void Haploid::transpose_mutate(Haploid& other, URNG& rng) {
    auto copying_transposons = this->transpose(rng);
    {
        auto tmp = other.transpose(rng);
        copying_transposons.insert(copying_transposons.end(),
            std::make_move_iterator(tmp.begin()),
            std::make_move_iterator(tmp.end()));
    }
    constexpr uint_fast32_t tolerance = 100;
    for (auto& p: copying_transposons) {
        for (uint_fast32_t i=0; i<tolerance; ++i) {
            auto target_haploid = this;
            if (rng.canonical() < 0.5) {
                target_haploid = &other;
            }
            auto& target_site = target_haploid->sites_[UNIFORM_SITES_(rng)];
            if (target_site) continue;
            target_site = std::move(p);
            break;
        }
    }
    this->mutate(rng);
    other.mutate(rng);
}

void Haploid::recombine(Haploid& other, URNG& rng) {
    bool flg = false;
    uint_fast32_t prev = 0U;
    for (auto this_it = this->sites_.begin(), other_it = other.sites_.begin();
         this_it != this->sites_.end() || other_it != other.sites_.end();
        ) {
        const uint_fast32_t this_pos = (this_it == this->end())?
            std::numeric_limits<uint_fast32_t>::max():
            this_it->first;
        const uint_fast32_t other_pos = (other_it == other.end())?
            std::numeric_limits<uint_fast32_t>::max():
            other_it->first;
        const uint_fast32_t here = std::min(this_pos, other_pos);
        std::poisson_distribution<uint_fast32_t> poisson((here - prev) * RECOMBINATION_RATE_);
        flg ^= (poisson(rng) % 2U);
        prev = here;
        if (this_pos < other_pos) {
            if (flg) {
                other.sites_.emplace_hint(other_it, std::move(*this_it));
                this_it = this->sites_.erase(this_it);
            } else {
                ++this_it;
            }
        } else {
            if (this_pos == other_pos) {
                if (flg) {
                  this_it->second.swap(other_it->second);
                }
                ++this_it;
                ++other_it;
            } else {
                if (flg) {
                    this->sites_.emplace_hint(this_it, std::move(*other_it));
                    other_it = other.sites_.erase(other_it);
                } else {
                    ++other_it;
                }
            }
        }
    }
}

void Haploid::mutate(URNG& rng) {
    using cnt_t = decltype(NUM_MUTATIONS_DIST_)::result_type;
    for (auto& p: sites_) {
        const cnt_t num_mutations = NUM_MUTATIONS_DIST_(rng);
        const bool is_deactivating = rng.canonical() < INDEL_RATE_;
        if (num_mutations > 0 || is_deactivating) {
            p.second = std::make_shared<Transposon>(*p.second);
        }
        for (cnt_t i=0; i<num_mutations; ++i) {
            p.second->mutate(rng);
        }
        if (is_deactivating) {
            p.second->indel();
        }
    }
}

void Haploid::evaluate_sites() {
    copy_number_ = sites_.size();
    prod_1_zs_ = 1.0;
    for (const auto& p: sites_) {
        prod_1_zs_ *= (1.0 - SELECTION_COEFS_GP_[p.first]);
    }
}

double Haploid::fitness(const Haploid& other) const {
    const double s_cn = XI_ * std::pow(copy_number_ + other.copy_number_, TAU_);
    return std::max(prod_1_zs_ * other.prod_1_zs_ * (1.0 - s_cn), 0.0);
}

std::vector<std::string> Haploid::summarize() const {
    // "site:indel:nonsynonymous:synonymous:activity"
    std::vector<std::string> v;
    v.reserve(copy_number_);
    for (const auto& p: sites_) {
        std::ostringstream oss;
        p.second->write_summary(oss << p.first << ":");
        v.push_back(oss.str());
    }
    return v;
}

std::map<double, uint_fast32_t> Haploid::count_activity() const {
    std::map<double, uint_fast32_t> counter;
    for (const auto& p: sites_) {
        ++counter[p.second->activity()];
    }
    return counter;
}

std::ostream& Haploid::write_positions(std::ostream& ost) const {
    return ost << wtl::str_join(sites_, ",", wtl::make_oss(), [](const auto& p){return p.first;});
}

std::ostream& Haploid::write_fasta(std::ostream& ost) const {
    for (const auto& p: sites_) {
        p.second->write_fasta(ost << ">" << this);
    }
    return ost;
}

std::ostream& operator<<(std::ostream& ost, const Haploid& x) {
    return ost << x.summarize();
}

void Haploid::test() {HERE;
    Haploid x = Haploid::copy_founder();
    std::cout << x << std::endl;
    x.write_fasta(std::cout);
    test_selection_coefs_gp();
    test_selection_coefs_cn();
    test_recombination();
}

void Haploid::test_selection_coefs_gp() {HERE;
    std::ofstream("tek-selection_coefs_gp.tsv")
        << "s_gp\n"
        << wtl::str_join(wtl::sample(SELECTION_COEFS_GP_, 2000, wtl::sfmt()), "\n");
    /*R
    read_tsv('tek-selection_coefs_gp.tsv') %>% {
      ggplot(., aes(s_gp))+
      geom_histogram(bins=30)+
      geom_vline(xintercept=mean(.$s_gp), colour='tomato')+
      theme_bw()
    } %>% {ggsave('selection_coefs_gp.pdf', ., width=4, height=4)}
    */
}

void Haploid::test_selection_coefs_cn() {HERE;
    std::ofstream ost("tek-selection_coefs_cn.tsv");
    ost << "xi\tcopy_number\ts_cn\n";
    const uint_fast32_t n = 2 * NUM_SITES;
    for (const double xi: {1e-5, 1e-4, 1e-3}) {
        for (uint_fast32_t i=0; i<n; ++i) {
            const double s_cn = xi * std::pow(i, TAU_);
            if (s_cn > 1.0) break;
            ost << xi << "\t" << i << "\t" << s_cn << "\n";
        }
    }
    /*R
    read_tsv('tek-selection_coefs_cn.tsv') %>%
    mutate(xi= sprintf('%.0e', xi)) %>% {
      ggplot(., aes(copy_number, s_cn, group=xi, colour=xi))+
      geom_line()+
      theme_bw()+theme(legend.position='top')
    } %>% {ggsave('selection_coefs_cn.pdf', ., width=4, height=4)}
    */
}

void Haploid::test_recombination() {HERE;
    Haploid zero;
    Haploid one;
    for (uint_fast32_t x=0U; x<60U; ++x) {
        one.sites_[x] = ORIGINAL_TE_;
    }
    zero.recombine(one, wtl::sfmt());
    one.write_positions(std::cerr) << std::endl;
    zero.write_positions(std::cerr) << std::endl;
}

} // namespace tek
