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
double Haploid::EXCISION_RATE_ = 1e-5;
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

/*! @ingroup params

    Command line option | Symbol        | Variable
    ------------------- | ------------- | -------------------------
    `--xi`              | \f$\xi\f$     | Haploid::XI_
    `--nu`              | \f$\nu\f$     | Haploid::EXCISION_RATE_
    `--lambda`          | \f$\lambda\f$ | Haploid::MEAN_SELECTION_COEF_
*/
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
    return founder;
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

Haploid Haploid::gametogenesis(const Haploid& other, URNG& rng) const {
    constexpr uint_fast32_t uint_max = std::numeric_limits<uint_fast32_t>::max();
    Haploid gamete(*this);
    bool flg = (rng.canonical() < 0.5);
    auto gamete_it = gamete.sites_.begin();
    auto gamete_end = gamete.sites_.end();
    auto other_it = other.sites_.cbegin();
    auto other_end = other.sites_.cend();
    uint_fast32_t gamete_pos = (gamete_it != gamete_end) ? gamete_it->first : uint_max;
    uint_fast32_t other_pos = (other_it != other_end) ? other_it->first : uint_max;
    uint_fast32_t here = 0U;
    uint_fast32_t prev = 0U;
    while ((here = std::min(gamete_pos, other_pos)) < uint_max) {
        std::poisson_distribution<uint_fast32_t> poisson((here - prev) * RECOMBINATION_RATE_);
        flg ^= (poisson(rng) % 2U);
        prev = here;
        if (gamete_pos < other_pos) {
            if (flg) {
                gamete_it = gamete.sites_.erase(gamete_it);
            } else {
                ++gamete_it;
            }
            gamete_pos = (gamete_it != gamete_end) ? gamete_it->first : uint_max;
        } else if (gamete_pos == other_pos) {
            if (flg) {
                gamete_it->second = other_it->second;
            }
            ++gamete_it;
            ++other_it;
            gamete_pos = (gamete_it != gamete_end) ? gamete_it->first : uint_max;
            other_pos = (other_it != other_end) ? other_it->first : uint_max;
        } else {
            if (flg) {
                gamete.sites_.emplace_hint(gamete_it, *other_it);
            }
            ++other_it;
            other_pos = (other_it != other_end) ? other_it->first : uint_max;
        }
    }
    return gamete;
}

std::vector<std::shared_ptr<Transposon>> Haploid::transpose(URNG& rng) {
    std::vector<std::shared_ptr<Transposon>> copying_transposons;
    for (auto it=sites_.cbegin(); it!=sites_.cend();) {
        if (rng.canonical() < it->second->transposition_rate()) {
            copying_transposons.push_back(it->second);
        }
        if (rng.canonical() < EXCISION_RATE_) {
            it = sites_.erase(it);
        } else {
            ++it;
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

double Haploid::prod_1_zs() const {
    double product = 1.0;
    for (const auto& p: sites_) {
        product *= (1.0 - SELECTION_COEFS_GP_[p.first]);
    }
    return product;
}

double Haploid::fitness(const Haploid& other) const {
    std::map<uint_fast32_t, uint_fast32_t> counter;
    for (const auto& p: this->sites_) {
        ++counter[p.second->species()];
    }
    for (const auto& p: other.sites_) {
        ++counter[p.second->species()];
    }
    double prod_1_xi_n_tau = 1.0;
    for (const auto& p: counter) {
        prod_1_xi_n_tau *= (1.0 - XI_ * std::pow(p.second, TAU_));
    }
    return std::max(prod_1_zs() * other.prod_1_zs() * prod_1_xi_n_tau, 0.0);
}

std::vector<std::string> Haploid::summarize() const {
    // "site:species:indel:nonsynonymous:synonymous:activity"
    std::vector<std::string> v;
    v.reserve(sites_.size());
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

//! shortcut of << Haploid::summarize()
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
    auto gamete = zero.gametogenesis(one, wtl::sfmt());
    gamete.write_positions(std::cerr) << std::endl;
}

} // namespace tek
