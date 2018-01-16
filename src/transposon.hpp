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
#include <unordered_map>
#include <random>
#include <mutex>

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

    Transposon(DNA&& non, DNA&& syn)
    : nonsynonymous_sites_(std::forward<DNA>(non)),
      synonymous_sites_(std::forward<DNA>(syn)) {}

    //! make one point mutation
    template <class URBG>
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
            speciate();
        }
    }

    //! set new #species_
    void speciate() {
        std::lock_guard<std::mutex> lock(MTX_);
        species_ = NUM_SPECIES_;
        ++NUM_SPECIES_;
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

    //! Hamming distance
    size_t operator-(const Transposon& other) const {
        return wtl::count(nonsynonymous_sites() != other.nonsynonymous_sites()) +
               wtl::count(synonymous_sites() != other.synonymous_sites());
    }
    //! interaction coefficient between species
    /*! \f[
            f(d) = \frac {L - d} {L - d _c}
        \f]
    */
    double operator*(const Transposon& other) const {
        const static double over_LEN_MIN_D = 1.0 / (LENGTH - MIN_DISTANCE());
        const auto dist = (*this - other);
        if (dist < MIN_DISTANCE()) {
            return 1.0;
        } else {
            return (LENGTH - dist) * over_LEN_MIN_D;
        }
    }
    //! clear #INTERACTION_COEFS_
    static void INTERACTION_COEFS_clear() {INTERACTION_COEFS_.clear();}
    //! setter of #INTERACTION_COEFS_
    static void INTERACTION_COEFS_emplace(uint_fast32_t x, uint_fast32_t y, double coef) {
        INTERACTION_COEFS_.emplace((static_cast<uint_fast64_t>(x) << 32) | y, coef);
    }
    //! getter of #INTERACTION_COEFS_
    static double INTERACTION_COEFS_get(uint_fast32_t x, uint_fast32_t y) {
        return INTERACTION_COEFS_.at((static_cast<uint_fast64_t>(x) << 32) | y);
    }
    //! getter of #INTERACTION_COEFS_
    static std::unordered_map<uint_fast64_t, double> INTERACTION_COEFS() {return INTERACTION_COEFS_;}
    //! getter of #MIN_DISTANCE_
    static size_t MIN_DISTANCE() {return MIN_DISTANCE_;}
    //! getter of #nonsynonymous_sites_
    const DNA& nonsynonymous_sites() const {return nonsynonymous_sites_;}
    //! getter of #synonymous_sites_
    const DNA& synonymous_sites() const {return synonymous_sites_;}
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
    std::ostream& write_fasta(std::ostream&, const uint_fast32_t copy_number=0u) const;
    //! write sequence
    std::ostream& write_sequence(std::ostream&) const;
    //! calculate and write activity for the given alpha and beta
    static void write_activity(std::ostream&, const double alpha, const unsigned int beta);
    friend std::ostream& operator<<(std::ostream&, const Transposon&);

    //! set dependent variables once, and reset non-parameter variables
    static void initialize();
    //! options description for Transposon class
    static boost::program_options::options_description options_desc();

  private:

    //!
    /*! \f[
            a = \left(\frac {1 - d - \alpha} {1 - \alpha} \right)^\beta
        \f]
    */
    static double calc_activity(const size_t num_mutations);

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    //! @addtogroup params
    //! @{

    //! \f$\alpha\f$, intercept of sequence identity to make activity zero
    static double ALPHA_;
    //! \f$\beta\f$, exponent of activity curve
    static unsigned int BETA_;
    //! speciation rate per mutation
    static double SPECIATION_RATE_;
    //! distance required for speciation
    static size_t MIN_DISTANCE_;
    //! @} params
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! 1 - #ALPHA_
    static double THRESHOLD_;
    //! pre-calculated activity values
    static std::array<double, NUM_NONSYNONYMOUS_SITES> ACTIVITY_;
    //! number of species; incremented by speciation
    static uint_fast32_t NUM_SPECIES_;
    //! interaction coefficients between species
    static std::unordered_map<uint_fast64_t, double> INTERACTION_COEFS_;
    //! mutex for speciation
    static std::mutex MTX_;

    //! nonsynonymous sites
    DNA nonsynonymous_sites_;
    //! synonymous sites
    DNA synonymous_sites_;
    //! activity is zero if this is true
    bool has_indel_ = false;
    //! transposon species
    uint_fast32_t species_ = 0;
};

/*! @brief TransposonFamily class
*/
class TransposonFamily {
  public:
    TransposonFamily()
    : nonsynonymous_sites_(Transposon::NUM_NONSYNONYMOUS_SITES),
      synonymous_sites_(Transposon::NUM_SYNONYMOUS_SITES) {}

    TransposonFamily& operator+=(const Transposon& x) {
        nonsynonymous_sites_ += x.nonsynonymous_sites();
        synonymous_sites_ += x.synonymous_sites();
        ++size_;
        return *this;
    }

    Transposon majority() const {
        return Transposon(nonsynonymous_sites_.majority(), synonymous_sites_.majority());
    }

    uint_fast32_t size() const {return size_;}

  private:
    Homolog nonsynonymous_sites_;
    Homolog synonymous_sites_;
    uint_fast32_t size_ = 0u;
};

} // namespace tek

#endif /* TEK_TRANSPOSON_HPP_ */
