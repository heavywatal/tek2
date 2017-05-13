// -*- mode: c++; coding: utf-8 -*-
/*! @file transposon.cpp
    @brief Implementation of Transposon class
*/
#include "transposon.hpp"

#include <iostream>

#include <wtl/debug.hpp>
#include <wtl/math.hpp>
#include <wtl/prandom.hpp>

namespace tek {

double Transposon::ALPHA_ = 0.8;
double Transposon::THRESHOLD_ = 1.0 - Transposon::ALPHA_;
unsigned int Transposon::BETA_ = 24U;

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
}

void Transposon::mutate() {
    static std::uniform_int_distribution<size_t> POS_DIST(0U, LENGTH - 1U);
    size_t pos = POS_DIST(wtl::sfmt());
    if (pos >= NUM_NONSYNONYMOUS_SITES_) {
        synonymous_sites_.flip(pos -= NUM_NONSYNONYMOUS_SITES_);
    } else {
        nonsynonymous_sites_.flip(pos);
    }
}

double Transposon::activity() const {
    if (has_indel_) return 0.0;
    const double diff = nonsynonymous_sites_.count() * OVER_NONSYNONYMOUS_SITES_;
    if (diff >= THRESHOLD_) return 0.0;
    return wtl::pow(1.0 - diff / THRESHOLD_, BETA_);
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

void Transposon::unit_test() {HERE;
    Transposon te;
    te.mutate();
    te.mutate();
    std::cout << te << std::endl;
    te.write_summary(std::cout) << std::endl;
    std::cout << te.activity() << std::endl;
    te.indel();
    te.write_summary(std::cout) << std::endl;
    std::cout << te.activity() << std::endl;
}

} // namespace tek
