/*! @file transposon.hpp
    @brief Interface of Transposon class
*/
#pragma once
#ifndef TEK_TRANSPOSON_HPP_
#define TEK_TRANSPOSON_HPP_

#include "dna.hpp"

#include <iosfwd>
#include <array>
#include <unordered_map>
#include <random>
#include <atomic>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

//! \f$L\f$, sequence length of TE (bp)
constexpr uint_fast32_t LENGTH = 300u;

//! @brief Parameters for Transposon class
/*! @ingroup params
*/
struct TransposonParams {
    //! \f$\alpha\f$, intercept of sequence identity to make activity zero
    double ALPHA = 0.7;
    //! \f$\beta\f$, exponent of activity curve
    unsigned int BETA = 6u;
    //! speciation rate per mutation
    double SPECIATION_RATE = 0.0;
    //! threshold distance required for speciation
    uint_fast32_t LOWER_THRESHOLD = LENGTH;
    //! threshold distance that interaction between sepecies becomes zero
    uint_fast32_t UPPER_THRESHOLD = LENGTH;
};

/*! @brief Transposon class
*/
class Transposon {
  public:
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    //! @addtogroup params
    //! @{
    //! Alias
    using param_type = TransposonParams;
    //! number of synonymous sites
    static constexpr uint_fast32_t NUM_SYNONYMOUS_SITES = LENGTH / 3u;
    //! number of nonsynonymous sites
    static constexpr uint_fast32_t NUM_NONSYNONYMOUS_SITES = LENGTH - NUM_SYNONYMOUS_SITES;
    //! \f$u_0\f$, maximum transposition rate
    static constexpr double MAX_TRANSPOSITION_RATE = 0.01;
    //! @} params
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! resiprocal of synonymous sites
    static constexpr double OVER_SYNONYMOUS_SITES = 1.0 / NUM_SYNONYMOUS_SITES;
    //! resiprocal of nonsynonymous sites
    static constexpr double OVER_NONSYNONYMOUS_SITES = 1.0 / NUM_NONSYNONYMOUS_SITES;

    //! default constructor
    Transposon() noexcept
    : nonsynonymous_sites_(NUM_NONSYNONYMOUS_SITES),
      synonymous_sites_(LENGTH - NUM_NONSYNONYMOUS_SITES) {}

    //! constructor
    Transposon(DNA&& non, DNA&& syn) noexcept
    : nonsynonymous_sites_(std::move(non)),
      synonymous_sites_(std::move(syn)) {}

    //! copy constructor
    Transposon(const Transposon&) = default;
    //! move constructor
    Transposon(Transposon&&) noexcept = default;

    //! make one point mutation
    template <class URBG>
    void mutate(URBG& engine) noexcept {
        thread_local std::uniform_int_distribution<uint_fast32_t> UNIF_LEN(0u, LENGTH - 1u);
        thread_local std::bernoulli_distribution BERN_SPECIATION(param().SPECIATION_RATE);
        auto pos = UNIF_LEN(engine);
        if (pos >= NUM_NONSYNONYMOUS_SITES) {
            synonymous_sites_.flip(pos -= NUM_NONSYNONYMOUS_SITES, engine);
        } else {
            nonsynonymous_sites_.flip(pos, engine);
        }
        if (param().SPECIATION_RATE > 0.0 && BERN_SPECIATION(engine)) {
            speciate();
        }
    }

    //! modify #species_ and #NUM_SPECIES_
    void speciate() noexcept {
        species_ = NUM_SPECIES_++;
    }

    //! set #has_indel_
    void indel() noexcept {has_indel_ = true;}

    //! \f$a_i\f$; count nonsynonymous mutations and return the pre-calculated #ACTIVITY_
    double activity() const noexcept {
        if (has_indel_) return 0.0;
        return ACTIVITY_[nonsynonymous_sites_.count()];
    }

    //! \f$u_i = u_0 \times a_i\f$
    double transposition_rate() const noexcept {
        return MAX_TRANSPOSITION_RATE * activity();
    }

    //! Hamming distance
    uint_fast32_t operator-(const Transposon& other) const noexcept {
        return (nonsynonymous_sites() - other.nonsynonymous_sites()) +
               (synonymous_sites() - other.synonymous_sites());
    }
    //! interaction coefficient between species
    /*! \f[
            I(d) = \begin{cases}
            1                            &         d < d_l \\
            \frac {d_u - d}{d_u - d_l}   & d_l \le d < d_u \\
            0                            & d_u \le d
            \end{cases}
        \f]
    */
    double operator*(const Transposon& other) const noexcept {
        const static double over_x = 1.0 / (param().UPPER_THRESHOLD - param().LOWER_THRESHOLD);
        const auto distance = (*this - other);
        if (distance < param().LOWER_THRESHOLD) {
            return 1.0;
        } else if (distance < param().UPPER_THRESHOLD) {
            return (param().UPPER_THRESHOLD - distance) * over_x;
        } else {
            return 0.0;
        }
    }
    //! check if distance is large enough for speciation
    bool is_far_enough_from(const Transposon& other) const noexcept {
        return (*this - other) >= param().LOWER_THRESHOLD;
    }
    //! check if speciation is allowed under the condition
    static bool can_speciate() noexcept {
        return param().LOWER_THRESHOLD < LENGTH;
    }
    //! clear #INTERACTION_COEFS_
    static void INTERACTION_COEFS_clear() noexcept {INTERACTION_COEFS_.clear();}
    //! setter of #INTERACTION_COEFS_
    static void INTERACTION_COEFS_emplace(uint_fast32_t x, uint_fast32_t y, double coef) noexcept {
        INTERACTION_COEFS_.emplace((static_cast<uint_fast64_t>(x) << 32) | y, coef);
    }
    //! getter of #INTERACTION_COEFS_
    static double INTERACTION_COEFS_get(uint_fast32_t x, uint_fast32_t y) {
        return INTERACTION_COEFS_.at((static_cast<uint_fast64_t>(x) << 32) | y);
    }
    //! getter of #INTERACTION_COEFS_
    static std::unordered_map<uint_fast64_t, double> INTERACTION_COEFS() noexcept {return INTERACTION_COEFS_;}
    //! getter of #nonsynonymous_sites_
    const DNA& nonsynonymous_sites() const noexcept {return nonsynonymous_sites_;}
    //! getter of #synonymous_sites_
    const DNA& synonymous_sites() const noexcept {return synonymous_sites_;}
    //! getter of #has_indel_
    bool has_indel() const noexcept {return has_indel_;}
    //! getter of #species_
    uint_fast32_t species() const noexcept {return species_;}
    //! nonsynonymous substitution per nonsynonymous site
    double dn() const noexcept {return nonsynonymous_sites_.count() * OVER_NONSYNONYMOUS_SITES;}
    //! synonymous substitution per synonymous site
    double ds() const noexcept {return synonymous_sites_.count() * OVER_SYNONYMOUS_SITES;}

    //! write summary
    std::ostream& write_summary(std::ostream&) const;
    //! write sequqnce with header in FASTA format
    std::ostream& write_fasta(std::ostream&) const;
    //! write metadata for FASTA header
    std::ostream& write_metadata(std::ostream&) const;
    //! write sequence
    std::ostream& write_sequence(std::ostream&) const;
    //! calculate and write activity for the given alpha and beta
    static void write_activity(std::ostream&, double alpha, unsigned int beta);
    friend std::ostream& operator<<(std::ostream&, const Transposon&);

    //! Set #PARAM_
    static void param(const param_type& p);
    //! Get #PARAM_
    static const param_type& param() {return PARAM_;}
    //! Set #PARAM_ with default values;
    static void initialize() {
        TransposonParams p;
        param(p);
    }

  private:
    //! Parameters shared among instances
    static param_type PARAM_;

    //!
    /*! \f[
            a = \left(\frac {1 - d - \alpha} {1 - \alpha} \right)^\beta
        \f]
    */
    static double calc_activity(uint_fast32_t num_mutations);

    //! 1 - #ALPHA_
    static double THRESHOLD_;
    //! pre-calculated activity values
    static std::array<double, NUM_NONSYNONYMOUS_SITES> ACTIVITY_;
    //! number of species; incremented by speciation
    static std::atomic_uint_fast32_t NUM_SPECIES_;
    //! interaction coefficients between species
    static std::unordered_map<uint_fast64_t, double> INTERACTION_COEFS_;

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
    TransposonFamily() noexcept
    : nonsynonymous_sites_(Transposon::NUM_NONSYNONYMOUS_SITES),
      synonymous_sites_(Transposon::NUM_SYNONYMOUS_SITES) {}

    TransposonFamily& operator+=(const Transposon& x) noexcept {
        nonsynonymous_sites_ += x.nonsynonymous_sites();
        synonymous_sites_ += x.synonymous_sites();
        ++size_;
        return *this;
    }

    Transposon majority() const noexcept {
        return Transposon(nonsynonymous_sites_.majority(), synonymous_sites_.majority());
    }

    uint_fast32_t size() const noexcept {return size_;}

  private:
    Homolog nonsynonymous_sites_;
    Homolog synonymous_sites_;
    uint_fast32_t size_ = 0u;
};

} // namespace tek

#endif /* TEK_TRANSPOSON_HPP_ */
