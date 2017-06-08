// -*- mode: c++; coding: utf-8 -*-
/*! @file transposon.cpp
    @brief Implementation of Transposon class
*/
#include "transposon.hpp"

#include <wtl/debug.hpp>
#include <wtl/prandom.hpp>

#include <cmath>
#include <iostream>
#include <fstream>

namespace tek {

double Transposon::ALPHA_ = 0.8;
double Transposon::THRESHOLD_ = 1.0 - Transposon::ALPHA_;
unsigned int Transposon::BETA_ = 24U;
std::array<double, Transposon::NUM_NONSYNONYMOUS_SITES> Transposon::ACTIVITY_;

namespace po = boost::program_options;

po::options_description Transposon::options_desc() {HERE;
    po::options_description description("Transposon");
    description.add_options()
      ("alpha", po::value(&ALPHA_)->default_value(ALPHA_))
      ("beta", po::value(&BETA_)->default_value(BETA_));
    return description;
}

void Transposon::set_parameters() {HERE;
    THRESHOLD_ = 1.0 - ALPHA_;
    for (size_t i=0; i<NUM_NONSYNONYMOUS_SITES; ++i) {
        ACTIVITY_[i] = calc_activity(i);
    }
}

double Transposon::calc_activity(const size_t n) {
    const double diff = n * OVER_NONSYNONYMOUS_SITES;
    if (diff >= THRESHOLD_) return 0.0;
    return std::pow(1.0 - diff / THRESHOLD_, BETA_);
}

void Transposon::mutate() {
    static std::uniform_int_distribution<size_t> POS_DIST(0U, LENGTH - 1U);
    size_t pos = POS_DIST(wtl::sfmt());
    if (pos >= NUM_NONSYNONYMOUS_SITES) {
        synonymous_sites_.flip(pos -= NUM_NONSYNONYMOUS_SITES);
    } else {
        nonsynonymous_sites_.flip(pos);
    }
}

std::ostream& Transposon::write_summary(std::ostream& ost) const {
    return ost << "(" << has_indel_ << ":"
               << nonsynonymous_sites_.count() << ":"
               << synonymous_sites_.count() << ")";
}

std::ostream& Transposon::write(std::ostream& ost) const {
    return ost << nonsynonymous_sites_ << synonymous_sites_;
}

std::ostream& operator<<(std::ostream& ost, const Transposon& x) {
    return x.write(ost);
}

void Transposon::test() {HERE;
    Transposon te;
    te.mutate();
    te.mutate();
    std::cout << te << std::endl;
    te.write_summary(std::cout) << std::endl;
    std::cout << te.activity() << std::endl;
    te.indel();
    te.write_summary(std::cout) << std::endl;
    std::cout << te.activity() << std::endl;
    test_activity();
}

void Transposon::test_activity() {HERE;
    std::ofstream ost("activity.tsv");
    ost << "alpha\tbeta\tidentity\tactivity\n";
    test_activity(ost, 0.70,  6);
    test_activity(ost, 0.75, 12);
    test_activity(ost, 0.80, 24);
    test_activity(ost, 0.85, 48);
    /* R
    read_tsv('activity.tsv') %>%
    {ggplot(., aes(identity, activity, group=alpha, colour=alpha)) + geom_line()} %>%
    {ggsave('fig1.pdf', ., width=5, height=3)}
    */
}

void Transposon::test_activity(std::ostream& ost, const double alpha, const unsigned int beta) {
    ALPHA_ = alpha;
    BETA_ = beta;
    set_parameters();
    for (size_t i=0; i<NUM_NONSYNONYMOUS_SITES; ++i) {
        double identity = (NUM_NONSYNONYMOUS_SITES - i) * OVER_NONSYNONYMOUS_SITES;
        if (identity < 0.7) continue;
        ost << alpha << "\t" << beta << "\t"
            << identity << "\t" << ACTIVITY_[i] << "\n";
    }
}

} // namespace tek
