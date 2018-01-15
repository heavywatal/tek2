/*! @file transposon.cpp
    @brief Implementation of Transposon class
*/
#include "transposon.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>

#include <cmath>
#include <iostream>

namespace tek {

double Transposon::ALPHA_ = 0.7;
unsigned int Transposon::BETA_ = 6u;
double Transposon::SPECIATION_RATE_ = 0.0;
size_t Transposon::MIN_DISTANCE_ = Transposon::LENGTH;
double Transposon::THRESHOLD_ = 0.0;
std::array<double, Transposon::NUM_NONSYNONYMOUS_SITES> Transposon::ACTIVITY_;
uint_fast32_t Transposon::NUM_SPECIES_ = 1u;
std::mutex Transposon::MTX_;

namespace po = boost::program_options;

/*! @ingroup params

    Command line option | Symbol        | Variable
    ------------------- | ------------- | -------------------------
    `-a,--alpha`        | \f$\alpha\f$  | Transposon::ALPHA_
    `-b,--beta`         | \f$\beta\f$   | Transposon::BETA_
    `--spec`            |               | Transposon::SPECIATION_RATE_
    `-d,--mindist`      |               | Transposon::MIN_DISTANCE_
*/
po::options_description Transposon::options_desc() {HERE;
    po::options_description description("Transposon");
    description.add_options()
      ("alpha,a", po::value(&ALPHA_)->default_value(ALPHA_))
      ("beta,b", po::value(&BETA_)->default_value(BETA_))
      ("spec", po::value(&SPECIATION_RATE_)->default_value(SPECIATION_RATE_))
      ("mindist,d", po::value(&MIN_DISTANCE_)->default_value(MIN_DISTANCE_));
    return description;
}

void Transposon::set_parameters() {HERE;
    THRESHOLD_ = 1.0 - ALPHA_;
    for (size_t i=0u; i<NUM_NONSYNONYMOUS_SITES; ++i) {
        ACTIVITY_[i] = calc_activity(i);
    }
}

double Transposon::calc_activity(const size_t num_mutations) {
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
    if (copy_number > 0u) ost << " copy_number=" << copy_number;
    return write_sequence(ost << "\n") << "\n";
}

std::ostream& Transposon::write_sequence(std::ostream& ost) const {
    for (size_t in=0u, is=0u; in<NUM_NONSYNONYMOUS_SITES; ++in, ++is) {
        ost << nonsynonymous_sites_.at(in);
        ost << nonsynonymous_sites_.at(++in);
        ost <<    synonymous_sites_.at(is);
    }
    return ost;
}

//! shortcut for Transposon::write_summary()
std::ostream& operator<<(std::ostream& ost, const Transposon& x) {
    return x.write_summary(ost);
}

void Transposon::write_activity(std::ostream& ost, const double alpha, const unsigned int beta) {
    ALPHA_ = alpha;
    BETA_ = beta;
    set_parameters();
    for (size_t i=0u; i<NUM_NONSYNONYMOUS_SITES; ++i) {
        double identity = (NUM_NONSYNONYMOUS_SITES - i) * OVER_NONSYNONYMOUS_SITES;
        if (identity < 0.7) continue;
        ost << alpha << "\t" << beta << "\t"
            << identity << "\t" << ACTIVITY_[i] << "\n";
    }
}

} // namespace tek
