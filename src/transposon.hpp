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
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    //! @addtogroup params
    //! @{

    //! \f$L\f$, sequence length of TE (bp)
    static constexpr size_t LENGTH = 300u;
    //! number of synonymous sites
    static constexpr size_t NUM_SYNONYMOUS_SITES = LENGTH / 3u;
    //! number of nonsynonymous sites
    static constexpr size_t NUM_NONSYNONYMOUS_SITES = LENGTH - NUM_SYNONYMOUS_SITES;
    //! \f$u_0\f$, maximum transposition rate
    static constexpr double MAX_TRANSPOSITION_RATE = 0.01;
    //! @} params
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! resiprocal of synonymous sites
    static constexpr double OVER_SYNONYMOUS_SITES = 1.0 / NUM_SYNONYMOUS_SITES;
    //! resiprocal of nonsynonymous sites
    static constexpr double OVER_NONSYNONYMOUS_SITES = 1.0 / NUM_NONSYNONYMOUS_SITES;

    //! default constructor
    Transposon()
    : nonsynonymous_sites_(NUM_NONSYNONYMOUS_SITES),
      synonymous_sites_(LENGTH - NUM_NONSYNONYMOUS_SITES) {}

    //! make one point mutation
    template <class URBG> inline
    void mutate(URBG& generator) {
        thread_local std::uniform_int_distribution<size_t> UNIF_LEN(0u, LENGTH - 1u);
        thread_local std::bernoulli_distribution BERN_SPECIATION(SPECIATION_RATE_);
        size_t pos = UNIF_LEN(generator);
        if (pos >= NUM_NONSYNONYMOUS_SITES) {
            synonymous_sites_.flip(pos -= NUM_NONSYNONYMOUS_SITES, generator);
        } else {
            nonsynonymous_sites_.flip(pos, generator);
        }
        if (SPECIATION_RATE_ > 0.0 && BERN_SPECIATION(generator)) {
            species_ = generator();
        }
    }

    //! set #has_indel_
    void indel() {has_indel_ = true;}
    //! \f$a_i\f$; count nonsynonymous mutations and return the pre-calculated #ACTIVITY_
    double activity() const {
        if (has_indel_) return 0.0;
        return ACTIVITY_[nonsynonymous_sites_.count()];
    }
    //! \f$u_i = u_0 \times a_i\f$
    double transposition_rate() const {
        return MAX_TRANSPOSITION_RATE * activity();
    }

    //! getter of #has_indel_
    bool has_indel() const {return has_indel_;}
    //! getter of #species_
    int32_t species() const {return species_;}
    //! nonsynonymous substitution per nonsynonymous site
    double dn() const {return nonsynonymous_sites_.count() * OVER_NONSYNONYMOUS_SITES;}
    //! synonymous substitution per synonymous site
    double ds() const {return synonymous_sites_.count() * OVER_SYNONYMOUS_SITES;}

    //! write summary
    std::ostream& write_summary(std::ostream&) const;
    //! write sequqnce with header in FASTA format
    std::ostream& write_fasta(std::ostream&, const uint_fast32_t copy_number=0u) const;
    //! write sequence
    std::ostream& write_sequence(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Transposon&);

    //! set #THRESHOLD_ and #ACTIVITY_
    static void set_parameters();
    //! options description for Transposon class
    static boost::program_options::options_description options_desc();
    //! unit test
    static void test();

  private:

    //!
    /*! \f[
            a = \left(\frac {1 - d - \alpha} {1 - \alpha} \right)^\beta
        \f]
    */
    static double calc_activity(const size_t num_mutations);
    //! Figure 1
    static void test_activity();
    //! implementation of test_activity()
    static void test_activity(std::ostream&, const double alpha, const unsigned int beta);

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    //! @addtogroup params
    //! @{

    //! \f$\alpha\f$, intercept of sequence identity to make activity zero
    static double ALPHA_;
    //! \f$\beta\f$, exponent of activity curve
    static unsigned int BETA_;
    //! speciation rate per mutation
    static double SPECIATION_RATE_;
    //! @} params
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! 1 - #ALPHA_
    static double THRESHOLD_;
    //! pre-calculated activity values
    static std::array<double, NUM_NONSYNONYMOUS_SITES> ACTIVITY_;

    //! nonsynonymous sites
    DNA nonsynonymous_sites_;
    //! synonymous sites
    DNA synonymous_sites_;
    //! activity is zero if this is true
    bool has_indel_ = false;
    //! transposon species
    int32_t species_ = 0;
};

} // namespace tek

#endif /* TEK_TRANSPOSON_HPP_ */
