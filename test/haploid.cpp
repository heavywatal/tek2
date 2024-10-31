#include "haploid.hpp"

#include <sfmt.hpp>
#include <wtl/iostr.hpp>

#include <random>
#include <iostream>
#include <fstream>

inline void selection_coefs_gp() {
    tek::Haploid::insert_coefs_gp(2000u);
    std::ofstream ofs("tek-selection_coefs_gp.tsv");
    ofs.exceptions(std::ios_base::failbit | std::ios_base::badbit);
    ofs << "s_gp\n";
    for (const auto& p: tek::Haploid::SELECTION_COEFS_GP()) {
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

inline void selection_coefs_cn() {
    std::ofstream ofs("tek-selection_coefs_cn.tsv");
    ofs.exceptions(std::ios_base::failbit | std::ios_base::badbit);
    ofs << "xi\tcopy_number\ts_cn\n";
    constexpr uint_fast32_t n = 10'000u;
    constexpr double tau = 1.5;
    for (const double xi: {1e-5, 1e-4, 1e-3}) {
        for (uint_fast32_t i=0u; i<n; ++i) {
            const double s_cn = xi * std::pow(i, tau);
            if (s_cn > 1.0) break;
            ofs << xi << "\t" << i << "\t" << s_cn << "\n";
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

template <class Iter> inline
std::vector<tek::Haploid::position_t> positions(Iter begin, Iter end, size_t n=0) {
    std::vector<tek::Haploid::position_t> x;
    x.reserve(n);
    while (begin != end) {
        x.push_back(begin->first);
        ++begin;
    }
    return x;
}

inline void recombination() {
    constexpr size_t n = 60u;
    tek::Haploid::URBG engine(std::random_device{}());
    tek::Haploid zero;
    tek::Haploid one(n, engine);
    auto gamete = zero.gametogenesis(one, engine);
    const auto one_pos = positions(one.begin(), one.end(), n);
    const auto gam_pos = positions(gamete.begin(), gamete.end(), n);
    std::cout << one_pos << std::endl;
    std::cout << gam_pos << std::endl;
    auto it = gam_pos.begin();
    for (const auto x: one_pos) {
        if (it != gam_pos.end() && x == *it) {
          std::cout << 1;
          ++it;
        } else {
          std::cout << 0;
        }
    }
    std::cout << "\n";
}

int main() {
    tek::Haploid::initialize(500u, 0.01, 20000);
    tek::Haploid x = tek::Haploid::copy_founder();
    std::cout << x << std::endl;
    x.write_fasta(std::cout);
    selection_coefs_gp();
    selection_coefs_cn();
    recombination();
    return 0;
}
