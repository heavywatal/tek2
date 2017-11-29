/*! @file haploid.cpp
    @brief Implementation of Haploid class
*/
#include "haploid.hpp"
#include "transposon.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/random.hpp>

#include <cmath>
#include <iostream>
#include <numeric>

namespace tek {

double Haploid::XI_ = 1e-4;
double Haploid::EXCISION_RATE_ = 1e-5;
double Haploid::MEAN_SELECTION_COEF_ = 1e-4;

double Haploid::MUTATION_RATE_ = 0.0;
double Haploid::RECOMBINATION_RATE_ = 0.0;
double Haploid::INDEL_RATE_ = 0.0;
std::unordered_map<Haploid::position_t, double> Haploid::SELECTION_COEFS_GP_;
std::shared_ptr<Transposon> Haploid::ORIGINAL_TE_ = std::make_shared<Transposon>();
std::shared_timed_mutex Haploid::MTX_;

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

void Haploid::set_parameters(const size_t popsize, const double theta, const double rho) {HERE;
    const double four_n = 4.0 * popsize;
    MUTATION_RATE_ = Transposon::LENGTH * theta / four_n;
    INDEL_RATE_ = MUTATION_RATE_ * INDEL_RATIO_;
    RECOMBINATION_RATE_ = rho / four_n;
    DCERR("MUTATION_RATE_ = " << MUTATION_RATE_ << std::endl);
    DCERR("INDEL_RATE_ = " << INDEL_RATE_ << std::endl);
    DCERR("RECOMBINATION_RATE_ = " << RECOMBINATION_RATE_ << std::endl);
}

Haploid::position_t Haploid::new_position(URBG& rng) {
    thread_local std::exponential_distribution<double> EXPO_DIST(1.0 / MEAN_SELECTION_COEF_);
    thread_local std::bernoulli_distribution BERN_FUNCTIONAL(PROP_FUNCTIONAL_SITES_);
    auto coef = BERN_FUNCTIONAL(rng) ? EXPO_DIST(rng) : 0.0;
    position_t j = 0u;
    std::lock_guard<std::shared_timed_mutex> lock(MTX_);
    while (!SELECTION_COEFS_GP_.emplace(j = rng(), coef).second) {;}
    return j;
}

Haploid Haploid::copy_founder() {
    constexpr position_t pos = std::numeric_limits<position_t>::max() / 2u;
    SELECTION_COEFS_GP_.emplace(pos, 0.0);
    Haploid founder;
    founder.sites_[pos] = ORIGINAL_TE_;
    return founder;
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

Haploid Haploid::gametogenesis(const Haploid& other, URBG& rng) const {
    constexpr position_t max_pos = std::numeric_limits<position_t>::max();
    Haploid gamete(*this);
    bool flg = (wtl::generate_canonical(rng) < 0.5);
    auto gamete_it = gamete.sites_.begin();
    auto gamete_end = gamete.sites_.end();
    auto other_it = other.sites_.cbegin();
    auto other_end = other.sites_.cend();
    position_t gamete_pos = (gamete_it != gamete_end) ? gamete_it->first : max_pos;
    position_t other_pos = (other_it != other_end) ? other_it->first : max_pos;
    position_t here = 0u;
    position_t prev = 0u;
    while ((here = std::min(gamete_pos, other_pos)) < max_pos) {
        const double lambda = (here - prev) * RECOMBINATION_RATE_;
        if (lambda < 10.0) {
            std::poisson_distribution<position_t> poisson(lambda);
            if (poisson(rng) % 2u != 0u) {
                flg = !flg;
            }
        } else {
            if (wtl::generate_canonical(rng) < 0.5) {
                flg = !flg;
            }
        }
        prev = here;
        if (gamete_pos < other_pos) {
            if (flg) {
                gamete_it = gamete.sites_.erase(gamete_it);
            } else {
                ++gamete_it;
            }
            gamete_pos = (gamete_it != gamete_end) ? gamete_it->first : max_pos;
        } else if (gamete_pos == other_pos) {
            if (flg) {
                gamete_it->second = other_it->second;
            }
            ++gamete_it;
            ++other_it;
            gamete_pos = (gamete_it != gamete_end) ? gamete_it->first : max_pos;
            other_pos = (other_it != other_end) ? other_it->first : max_pos;
        } else {
            if (flg) {
                gamete.sites_.emplace_hint(gamete_it, *other_it);
            }
            ++other_it;
            other_pos = (other_it != other_end) ? other_it->first : max_pos;
        }
    }
    return gamete;
}

std::vector<std::shared_ptr<Transposon>> Haploid::transpose(URBG& rng) {
    std::vector<std::shared_ptr<Transposon>> copying_transposons;
    for (auto it=sites_.cbegin(); it!=sites_.cend();) {
        if (wtl::generate_canonical(rng) < it->second->transposition_rate()) {
            copying_transposons.push_back(it->second);
        }
        if (wtl::generate_canonical(rng) < EXCISION_RATE_) {
            it = sites_.erase(it);
        } else {
            ++it;
        }
    }
    return copying_transposons;
}

void Haploid::transpose_mutate(Haploid& other, URBG& rng) {
    auto copying_transposons = this->transpose(rng);
    {
        auto tmp = other.transpose(rng);
        copying_transposons.insert(copying_transposons.end(),
            std::make_move_iterator(tmp.begin()),
            std::make_move_iterator(tmp.end()));
    }
    for (auto& p: copying_transposons) {
        auto target_haploid = this;
        if (wtl::generate_canonical(rng) < 0.5) {
            target_haploid = &other;
        }
        target_haploid->sites_.emplace(new_position(rng), std::move(p));
    }
    this->mutate(rng);
    other.mutate(rng);
}

void Haploid::mutate(URBG& rng) {
    thread_local std::poisson_distribution<uint_fast32_t> POISSON_MUT(MUTATION_RATE_);
    thread_local std::bernoulli_distribution BERN_INDEL(INDEL_RATE_);
    for (auto& p: sites_) {
        const uint_fast32_t num_mutations = POISSON_MUT(rng);
        const bool is_deactivating = BERN_INDEL(rng);
        if (num_mutations > 0u || is_deactivating) {
            p.second = std::make_shared<Transposon>(*p.second);
        }
        for (uint_fast32_t i=0u; i<num_mutations; ++i) {
            p.second->mutate(rng);
        }
        if (is_deactivating) {
            p.second->indel();
        }
    }
}

double Haploid::prod_1_zs() const {
    double product = 1.0;
    std::shared_lock<std::shared_timed_mutex> lock(MTX_);
    for (const auto& p: sites_) {
        product *= (1.0 - SELECTION_COEFS_GP_.at(p.first));
    }
    return product;
}

double Haploid::fitness(const Haploid& other) const {
    std::map<int32_t, uint_fast32_t> counter;
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
    URBG rng(std::random_device{}());
    for (size_t i=SELECTION_COEFS_GP_.size(); i<2000u; ++i) {
        new_position(rng);
    }
    auto ofs = wtl::make_ofs("tek-selection_coefs_gp.tsv");
    ofs << "s_gp\n";
    for (const auto& p: SELECTION_COEFS_GP_) {
        ofs << p.second << "\n";
    }
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
    auto ost = wtl::make_ofs("tek-selection_coefs_cn.tsv");
    ost << "xi\tcopy_number\ts_cn\n";
    const uint_fast32_t n = 10'000u;
    for (const double xi: {1e-5, 1e-4, 1e-3}) {
        for (uint_fast32_t i=0u; i<n; ++i) {
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
    URBG rng(std::random_device{}());
    Haploid zero;
    Haploid one;
    for (uint_fast32_t x=0u; x<60u; ++x) {
        one.sites_[x] = ORIGINAL_TE_;
    }
    auto gamete = zero.gametogenesis(one, rng);
    gamete.write_positions(std::cerr) << std::endl;
}

} // namespace tek
