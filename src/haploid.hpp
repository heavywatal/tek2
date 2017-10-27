/*! @file haploid.hpp
    @brief Interface of Haploid class
*/
#pragma once
#ifndef TEK_HAPLOID_HPP_
#define TEK_HAPLOID_HPP_

#include <boost/program_options.hpp>

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>
#include <valarray>
#include <map>
#include <memory>
#include <random>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace wtl {
    class sfmt19937;
}

namespace tek {

class Transposon;

/*! @brief Haploid class
*/
class Haploid {
  public:
    //! random number generator class
    using URNG = wtl::sfmt19937;

    //! default constructor
    Haploid() = default;
    //! default copy constructor
    Haploid(const Haploid&) = default;
    //! default move constructor
    Haploid(Haploid&& other) = default;
    //! default move assignment operator
    Haploid& operator=(Haploid&&) = default;

    //! return a Haploid object after recombination
    Haploid gametogenesis(const Haploid& other, URNG& rng) const;
    //! mutation process within an individual
    void transpose_mutate(Haploid& other, URNG& rng);
    //! evaluate and return fitness
    /*! \f[\begin{split}
            s_{CN,k} &= \xi n_k ^\tau \\
            w_k &= w_{k,GP} (1 - s_{CN,k})
        \end{split}\f]
    */
    double fitness(const Haploid&) const;

    //! !sites_.empty()
    bool has_transposon() const {return !sites_.empty();};
    //! count activity
    std::map<double, uint_fast32_t> count_activity() const;

    //! return vector of Transposon summaries
    std::vector<std::string> summarize() const;
    //! write positions
    std::ostream& write_positions(std::ostream&) const;
    //! write sequence with address as name
    std::ostream& write_fasta(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Haploid&);

    //! shortcut of sites_.begin()
    auto begin() const {return sites_.begin();}
    //! shortcut of sites_.end()
    auto end() const {return sites_.end();}

    //! return a Haploid with an #ORIGINAL_TE_ on the same site
    static Haploid copy_founder();
    //! set static member variables
    static void set_parameters(const size_t popsize, const double theta, const double rho);
    //! options description for Haploid class
    static boost::program_options::options_description options_desc();
    //! unit test
    static void test();

  private:
    //! default copy assignment operator (private)
    Haploid& operator=(const Haploid&) = default;

    //! set \f$s_{GP}\f$ for all TE sites
    static void set_SELECTION_COEFS_GP();

    //! return TEs to be transposed
    std::vector<std::shared_ptr<Transposon>> transpose(URNG&);
    //! make point mutation, indel, and speciation
    void mutate(URNG&);
    //! calculate genome position component of fitness
    /*! \f[
            w_{k,GP} = \prod _j^T (1 - z_j s_{GP,j})
        \f]
    */
    double prod_1_zs() const;

    //! write to file
    static void test_selection_coefs_cn();
    //! write to file
    static void test_selection_coefs_gp();
    //! print to std::cerr
    static void test_recombination();

    //! @ingroup params
    //! \f$T\f$, number of TE sites in a haploid genome
    static constexpr uint_fast32_t NUM_SITES = 2'000'000u;
    //! @ingroup params
    //! \f$\phi\f$, relative rate of indels to point mutation
    static constexpr double INDEL_RATIO_ = 0.2;
    //! @ingroup params
    //! \f$\tau\f$, constant for the intensity of copy number selection
    static constexpr double TAU_ = 1.5;
    //! @ingroup params
    //! \f$p\f$, proportion of non-neutral sites
    static constexpr double PROP_FUNCTIONAL_SITES_ = 0.75;
    //! \f$\xi\f$, parameter for the intensity of copy number selection
    static double XI_;
    //! \f$\nu\f$, excision rate per generation per element
    static double EXCISION_RATE_;
    //! \f$\lambda\f$, mean selection coef against TEs on functional sites
    static double MEAN_SELECTION_COEF_;

    //! \f$c = \rho / 4N\f$, recombination rate per site
    static double RECOMBINATION_RATE_;
    //! \f$\phi\mu\f$, absolute indel rate
    static double INDEL_RATE_;
    //! pre-calculated coefficient of GP selection
    static std::valarray<double> SELECTION_COEFS_GP_;
    //! uniform distribution to get mutating site
    static std::uniform_int_distribution<uint_fast32_t> UNIFORM_SITES_;
    //! poisson distribution to get number of mutations
    static std::poisson_distribution<uint_fast32_t> NUM_MUTATIONS_DIST_;
    //! original TE with no mutation and complete activity
    static std::shared_ptr<Transposon> ORIGINAL_TE_;

    //! position => shptr to transposon
    std::map<uint_fast32_t, std::shared_ptr<Transposon>> sites_;
};

} // namespace tek

#endif /* TEK_HAPLOID_HPP_ */
