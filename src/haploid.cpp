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
std::uniform_int_distribution<size_t> Haploid::UNIFORM_SITES_(0, Haploid::NUM_SITES - 1U);
std::poisson_distribution<unsigned int> Haploid::NUM_MUTATIONS_DIST_(0.0);
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
    std::vector<size_t> sites(NUM_SITES);
    std::iota(sites.begin(), sites.end(), 0);
    const size_t n = NUM_SITES * PROP_FUNCTIONAL_SITES_;
    const auto functional_sites = wtl::sample(sites, n, wtl::sfmt());
    std::exponential_distribution<double> expo_dist(1.0 / MEAN_SELECTION_COEF_);
    for (const auto i: functional_sites) {
        SELECTION_COEFS_GP_[i] = expo_dist(wtl::sfmt());
    }
}

void Haploid::set_parameters(const size_t popsize, const double theta, const double rho) {HERE;
    const double four_n = 4.0 * popsize;
    RECOMBINATION_RATE_ = rho / four_n;
    std::cerr << "RECOMBINATION_RATE_ = " << RECOMBINATION_RATE_ << std::endl;
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
        if (!p) continue;
        if (rng.canonical() < p->transposition_rate()) {
            copying_transposons.push_back(p);
        }
        if (rng.canonical() < EXCISION_RATE_) {
            p.reset();
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
    constexpr unsigned int tolerance = 100;
    for (auto& p: copying_transposons) {
        for (unsigned int i=0; i<tolerance; ++i) {
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
    for (size_t j=1; j<NUM_SITES; ++j) {
        if (rng.canonical() < RECOMBINATION_RATE_) {
            flg = !flg;
        }
        if (flg) {
            sites_[j].swap(other.sites_[j]);
        }
    }
}

void Haploid::mutate(URNG& rng) {
    using cnt_t = decltype(NUM_MUTATIONS_DIST_)::result_type;
    for (auto& p: sites_) {
        if (!p) continue;
        const cnt_t num_mutations = NUM_MUTATIONS_DIST_(rng);
        const bool is_deactivating = rng.canonical() < INDEL_RATE_;
        if (num_mutations > 0 || is_deactivating) {
            p = std::make_shared<Transposon>(*p);
        }
        for (cnt_t i=0; i<num_mutations; ++i) {
            p->mutate(rng);
        }
        if (is_deactivating) {
            p->indel();
        }
    }
}

void Haploid::evaluate_sites() {
    copy_number_ = 0U;
    prod_1_zs_ = 1.0;
    for (size_t j=0; j<NUM_SITES; ++j) {
        if (sites_[j]) {
            ++copy_number_;
            prod_1_zs_ *= (1.0 - SELECTION_COEFS_GP_[j]);
        }
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
    for (size_t i=0; i<NUM_SITES; ++i) {
        if (sites_[i]) {
            std::ostringstream oss;
            sites_[i]->write_summary(oss << i << ":");
            v.push_back(oss.str());
        }
    }
    return v;
}

std::ostream& Haploid::write_binary(std::ostream& ost) const {
    for (const auto& p: sites_) {
        if (p) {ost << '1';}
        else   {ost << '0';}
    }
    return ost;
}

std::ostream& Haploid::write_fasta(std::ostream& ost) const {
    for (const auto& p: sites_) {
        if (p) p->write_fasta(ost << ">" << this);
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

void Haploid::test_selection_coefs_gp() {
    std::ofstream("tek-selection_coefs_gp.tsv")
        << "s_gp\n"
        << wtl::str_join(SELECTION_COEFS_GP_, "\n");
    /*R
    read_tsv('tek-selection_coefs_gp.tsv') %>% {
      ggplot(., aes(s_gp))+
      geom_histogram(bins=30)+
      geom_vline(xintercept=mean(.$s_gp), colour='tomato')+
      theme_bw()
    } %>% {ggsave('selection_coefs_gp.pdf', ., width=4, height=4)}
    */
}

void Haploid::test_selection_coefs_cn() {
    std::ofstream ost("tek-selection_coefs_cn.tsv");
    ost << "xi\tcopy_number\ts_cn\n";
    const size_t n = 2 * NUM_SITES;
    for (const double xi: {1e-5, 1e-4, 1e-3}) {
        for (size_t i=0; i<n; ++i) {
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

void Haploid::test_recombination() {
    Haploid zero;
    Haploid one;
    for (auto& p: one.sites_) {p = ORIGINAL_TE_;}
    zero.recombine(one, wtl::sfmt());
    std::ofstream ost("tek-recombination.txt");
    zero.write_binary(ost) << std::endl;
    /*R
    read_lines('tek-recombination.txt') %>%
    str_match_all(c('0+', '1+')) %>%
    flatten_chr() %>%
    {tibble(length=nchar(.))} %>% {
      ggplot(., aes(length))+
      geom_bar()+
      geom_vline(xintercept=mean(.$length), colour='tomato')+
      theme_bw()
    } %>% {ggsave('tract_length.pdf', ., width=4, height=4)}
    */
}

} // namespace tek
