/*! @file transposon.cpp
    @brief Implementation of Transposon class
*/
#include "transposon.hpp"

#include <wtl/debug.hpp>
#include <wtl/numeric.hpp>

#include <cmath>

namespace tek {

static_assert(std::is_nothrow_default_constructible<Transposon>{}, "");
static_assert(std::is_nothrow_move_constructible<Transposon>{}, "");

static_assert(std::is_nothrow_default_constructible<DNA<3>>{}, "");
static_assert(std::is_nothrow_move_constructible<DNA<3>>{}, "");

void Transposon::param(const param_type& p) {HERE;
    PARAM_ = p;
    static bool has_been_executed = false;
    NUM_SPECIES_.store(1u);
    if (has_been_executed) return;
    THRESHOLD_ = 1.0 - param().ALPHA;
    for (uint_fast32_t i=0u; i<NUM_NONSYNONYMOUS_SITES; ++i) {
        ACTIVITY_[i] = calc_activity(i);
    }
    has_been_executed = true;
}

double Transposon::calc_activity(uint_fast32_t num_mutations) {
    const double diff = num_mutations * OVER_NONSYNONYMOUS_SITES;
    if (diff >= THRESHOLD_) return 0.0;
    return std::pow(1.0 - diff / THRESHOLD_, param().BETA);
}

std::ostream& Transposon::write_summary(std::ostream& ost) const {
    return ost << species_ << ":"
               << has_indel_ << ":"
               << nonsynonymous_sites_.count() << ":"
               << synonymous_sites_.count() << ":"
               << activity();
}

std::ostream& Transposon::write_fasta(std::ostream& ost) const {
    write_metadata(ost << ">");
    write_sequence(ost << "\n");
    return ost << "\n";
}

std::ostream& Transposon::write_metadata(std::ostream& ost) const {
    return ost << "te=" << this
        << " species=" << species_ << " indel=" << has_indel_
        << " dn=" << dn() << " ds=" << ds()
        << " activity=" << activity();
}

std::ostream& Transposon::write_sequence(std::ostream& ost) const {
    for (uint_fast32_t in=0u, is=0u; in<NUM_NONSYNONYMOUS_SITES; ++in, ++is) {
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

void Transposon::write_activity(std::ostream& ost, const double alpha, const unsigned int beta) {
    for (uint_fast32_t i=0u; i<NUM_NONSYNONYMOUS_SITES; ++i) {
        double identity = (NUM_NONSYNONYMOUS_SITES - i) * OVER_NONSYNONYMOUS_SITES;
        if (identity < 0.7) continue;
        ost << alpha << "\t" << beta << "\t"
            << identity << "\t" << calc_activity(i) << "\n";
    }
}

} // namespace tek
