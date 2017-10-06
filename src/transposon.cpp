// -*- mode: c++; coding: utf-8 -*-
/*! @file transposon.cpp
    @brief Implementation of Transposon class
*/
#include "transposon.hpp"

#include <wtl/debug.hpp>

#include <cmath>
#include <iostream>
#include <fstream>

namespace tek {

double Transposon::ALPHA_ = 0.7;
unsigned int Transposon::BETA_ = 6U;
double Transposon::SPECIATION_RATE_ = 0.0;
double Transposon::THRESHOLD_ = 0.0;
std::array<double, Transposon::NUM_NONSYNONYMOUS_SITES> Transposon::ACTIVITY_;

namespace po = boost::program_options;

/*! @ingroup params

    Command line option | Symbol        | Variable
    ------------------- | ------------- | -------------------------
    `--alpha`           | \f$\alpha\f$  | Transposon::ALPHA_
    `--beta`            | \f$\beta\f$   | Transposon::BETA_
    `--spec`            |               | Transposon::SPECIATION_RATE_
*/
po::options_description Transposon::options_desc() {HERE;
    po::options_description description("Transposon");
    description.add_options()
      ("alpha", po::value(&ALPHA_)->default_value(ALPHA_))
      ("beta", po::value(&BETA_)->default_value(BETA_))
      ("spec", po::value(&SPECIATION_RATE_)->default_value(SPECIATION_RATE_));
    return description;
}

void Transposon::set_parameters() {HERE;
    THRESHOLD_ = 1.0 - ALPHA_;
    for (uint_fast32_t i=0; i<NUM_NONSYNONYMOUS_SITES; ++i) {
        ACTIVITY_[i] = calc_activity(i);
    }
}

double Transposon::calc_activity(const uint_fast32_t num_mutations) {
    const double diff = num_mutations * OVER_NONSYNONYMOUS_SITES;
    if (diff >= THRESHOLD_) return 0.0;
    return std::pow(1.0 - diff / THRESHOLD_, BETA_);
}

std::ostream& Transposon::write_summary(std::ostream& ost) const {
    return ost << species_ << ":"
               << has_indel_ << ":"
               << nonsynonymous_sites_.count() << ":"
               << synonymous_sites_.count() << ":"
               << activity();
}

std::ostream& Transposon::write_fasta(std::ostream& ost, const uint_fast32_t copy_number) const {
    ost << ">" << this
        << " species=" << species_ << " indel=" << has_indel_
        << " dn=" << dn() << " ds=" << ds()
        << " activity=" << activity();
    if (copy_number > 0U) ost << " copy_number=" << copy_number;
    return write_sequence(ost << "\n") << "\n";
}

std::ostream& Transposon::write_sequence(std::ostream& ost) const {
    for (uint_fast32_t in=0, is=0; in<NUM_NONSYNONYMOUS_SITES; ++in, ++is) {
        ost << nonsynonymous_sites_[in];
        ost << nonsynonymous_sites_[++in];
        ost <<    synonymous_sites_[is];
    }
    return ost;
}

//! shortcut for Transposon::write_summary()
std::ostream& operator<<(std::ostream& ost, const Transposon& x) {
    return x.write_summary(ost);
}

void Transposon::test() {HERE;
    Transposon te;
    std::mt19937 mt(std::random_device{}());
    te.mutate(mt);
    te.mutate(mt);
    te.write_fasta(std::cout);
    std::cout << te << std::endl;
    te.indel();
    std::cout << te << std::endl;
    test_activity();
}

void Transposon::test_activity() {HERE;
    std::ofstream ost("tek-activity_function.tsv");
    ost << "alpha\tbeta\tidentity\tactivity\n";
    test_activity(ost, 0.70,  6);
    test_activity(ost, 0.75, 12);
    test_activity(ost, 0.80, 24);
    test_activity(ost, 0.85, 48);
    /* R
    read_tsv('tek-activity_function.tsv') %>%
    {ggplot(., aes(identity, activity, group=alpha, colour=alpha)) + geom_line()} %>%
    {ggsave('activity_function.pdf', ., width=5, height=3)}
    */
}

void Transposon::test_activity(std::ostream& ost, const double alpha, const unsigned int beta) {
    ALPHA_ = alpha;
    BETA_ = beta;
    set_parameters();
    for (uint_fast32_t i=0; i<NUM_NONSYNONYMOUS_SITES; ++i) {
        double identity = (NUM_NONSYNONYMOUS_SITES - i) * OVER_NONSYNONYMOUS_SITES;
        if (identity < 0.7) continue;
        ost << alpha << "\t" << beta << "\t"
            << identity << "\t" << ACTIVITY_[i] << "\n";
    }
}

} // namespace tek
