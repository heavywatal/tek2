#include "haploid.hpp"

#include <wtl/iostr.hpp>

#include <random>
#include <iostream>

inline void selection_coefs_gp() {
    tek::Haploid::insert_coefs_gp(2000u);
    auto ofs = wtl::make_ofs("tek-selection_coefs_gp.tsv");
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
    auto ost = wtl::make_ofs("tek-selection_coefs_cn.tsv");
    ost << "xi\tcopy_number\ts_cn\n";
    constexpr uint_fast32_t n = 10'000u;
    constexpr double tau = 1.5;
    for (const double xi: {1e-5, 1e-4, 1e-3}) {
        for (uint_fast32_t i=0u; i<n; ++i) {
            const double s_cn = xi * std::pow(i, tau);
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

inline void recombination() {
    tek::Haploid::URBG engine(std::random_device{}());
    tek::Haploid zero;
    tek::Haploid one(60u);
    auto gamete = zero.gametogenesis(one, engine);
    gamete.write_positions(std::cerr) << std::endl;
}

int main() {
    tek::Haploid::initialize(500u, 0.01, 200);
    tek::Haploid x = tek::Haploid::copy_founder();
    std::cout << x << std::endl;
    x.write_fasta(std::cout);
    selection_coefs_gp();
    selection_coefs_cn();
    recombination();
    return 0;
}
