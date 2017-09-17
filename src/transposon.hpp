// -*- mode: c++; coding: utf-8 -*-
/*! @file transposon.hpp
    @brief Interface of Transposon class
*/
#pragma once
#ifndef TEK_TRANSPOSON_HPP_
#define TEK_TRANSPOSON_HPP_

#include "dna.hpp"

#include <boost/program_options.hpp>

#include <iosfwd>
#include <array>
#include <random>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

/*! @brief Transposon class
*/
class Transposon {
  public:
    //! number of base pairs
    static constexpr uint_fast32_t LENGTH = 300;
    //! number of synonymous sites
    static constexpr uint_fast32_t NUM_SYNONYMOUS_SITES = LENGTH / 3;
    //! number of nonsynonymous sites
    static constexpr uint_fast32_t NUM_NONSYNONYMOUS_SITES = LENGTH - NUM_SYNONYMOUS_SITES;
    //! resiprocal of synonymous sites
    static constexpr double OVER_SYNONYMOUS_SITES = 1.0 / NUM_SYNONYMOUS_SITES;
    //! resiprocal of nonsynonymous sites
    static constexpr double OVER_NONSYNONYMOUS_SITES = 1.0 / NUM_NONSYNONYMOUS_SITES;
    //! basal value of transposition rate
    static constexpr double MAX_TRANSPOSITION_RATE = 0.01;

    //! default constructor
    Transposon() = default;

    //! make one point mutation
    template <class URNG> inline
    void mutate(URNG& generator) {
        static std::uniform_int_distribution<uint_fast32_t> UNIF_LEN(0U, LENGTH - 1U);
        uint_fast32_t pos = UNIF_LEN(generator);
        if (pos >= NUM_NONSYNONYMOUS_SITES) {
            synonymous_sites_.flip(pos -= NUM_NONSYNONYMOUS_SITES, generator);
        } else {
            nonsynonymous_sites_.flip(pos, generator);
        }
    }

    //! modify #species_
    template <class URNG> inline
    void speciate(URNG& generator) {species_ = generator();}

    //! set #has_indel_
    void indel() {has_indel_ = true;}
    //! count nonsynonymous mutations and return one of #ACTIVITY_
    double activity() const {
        if (has_indel_) return 0.0;
        return ACTIVITY_[nonsynonymous_sites_.count()];
    }
    //! transposition rate
    double transposition_rate() const {
        return MAX_TRANSPOSITION_RATE * activity();
    }

    //! getter of #has_indel_
    bool has_indel() const {return has_indel_;}
    //! getter of #species_
    uint_fast32_t species() const {return species_;}
    //! nonsynonymous substitution per nonsynonymous site
    double dn() const {return nonsynonymous_sites_.count() * OVER_NONSYNONYMOUS_SITES;}
    //! synonymous substitution per synonymous site
    double ds() const {return synonymous_sites_.count() * OVER_SYNONYMOUS_SITES;}

    //! write summary
    std::ostream& write_summary(std::ostream&) const;
    //! write sequqnce with header in FASTA format
    std::ostream& write_fasta(std::ostream&, const uint_fast32_t copy_number=0) const;
    //! write sequence
    std::ostream& write_sequence(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Transposon&);

    //! set #THRESHOLD_ and #ACTIVITY_
    static void set_parameters();
    //! options description for optional arguments
    static boost::program_options::options_description options_desc();
    //! unit test
    static void test();

  private:

    //!
    /*! \f[
            a = \left(\frac {1 - d - \alpha} {1 - \alpha} \right)^\beta
        \f]
    */
    static double calc_activity(const uint_fast32_t num_mutations);
    //! Figure 1
    static void test_activity();
    //! implementation of test_activity()
    static void test_activity(std::ostream&, const double alpha, const unsigned int beta);

    //! intercept of sequence identity to make activity zero
    static double ALPHA_;
    //! 1 - #ALPHA_
    static double THRESHOLD_;
    //! exponent of activity curve
    static unsigned int BETA_;
    //! pre-calculated activity values
    static std::array<double, NUM_NONSYNONYMOUS_SITES> ACTIVITY_;

    //! nonsynonymous sites
    DNA<NUM_NONSYNONYMOUS_SITES> nonsynonymous_sites_;
    //! synonymous sites
    DNA<LENGTH - NUM_NONSYNONYMOUS_SITES> synonymous_sites_;
    //! activity is zero if this is true
    bool has_indel_ = false;
    //! transposon species
    uint_fast32_t species_ = 0;
};

} // namespace tek

#endif /* TEK_TRANSPOSON_HPP_ */
