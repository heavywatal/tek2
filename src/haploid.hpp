// -*- mode: c++; coding: utf-8 -*-
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
    using URNG = wtl::sfmt19937;

    //! default constructor
    Haploid() = default;
    //! default copy constructor
    Haploid(const Haploid&) = default;
    //! default move constructor
    Haploid(Haploid&& other) = default;
    //! default move assignment operator
    Haploid& operator=(Haploid&&) = default;

    Haploid gametogenesis(const Haploid& other, URNG& rng) const;
    void transpose_mutate(Haploid& other, URNG& rng);
    double fitness(const Haploid&) const;

    bool has_transposon() const {return !sites_.empty();};
    std::map<double, uint_fast32_t> count_activity() const;

    std::vector<std::string> summarize() const;
    std::ostream& write_positions(std::ostream&) const;
    std::ostream& write_fasta(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Haploid&);

    auto begin() const {return sites_.begin();}
    auto end() const {return sites_.end();}

    static Haploid copy_founder();
    static void set_parameters(const size_t popsize, const double theta, const double rho);
    static boost::program_options::options_description options_desc();
    static void test();

  private:
    Haploid& operator=(const Haploid&) = default;

    static void set_SELECTION_COEFS_GP();

    std::vector<std::shared_ptr<Transposon>> transpose(URNG&);
    void mutate(URNG&);
    double prod_1_zs() const;

    static void test_selection_coefs_cn();
    static void test_selection_coefs_gp();
    static void test_recombination();

    //! number of TE sites in a haploid genome \f$T\f$
    static constexpr uint_fast32_t NUM_SITES = 2'000'000;
    //! relative rate of indels to point mutation \f$\phi\f
    static constexpr double INDEL_RATIO_ = 0.2;
    //! constant for the intensity of copy number selection \f$\tau\f$
    static constexpr double TAU_ = 1.5;
    //! proportion of non-neutral sites
    static constexpr double PROP_FUNCTIONAL_SITES_ = 0.75;
    //! parameter for the intensity of copy number selection \f$\xi\f$
    static double XI_;
    //! excision rate per generation per element \f$\nu\f$
    static double EXCISION_RATE_;
    //! mean selection coef against TEs on functional sites \f$\lambda\f
    static double MEAN_SELECTION_COEF_;
    static double SPECIATION_RATIO_;

    //! recombination rate per site \f$c = \rho / 4N\f$
    static double RECOMBINATION_RATE_;
    static double INDEL_RATE_;
    static double SPECIATION_RATE_;
    static std::valarray<double> SELECTION_COEFS_GP_;
    static std::uniform_int_distribution<uint_fast32_t> UNIFORM_SITES_;
    static std::poisson_distribution<uint_fast32_t> NUM_MUTATIONS_DIST_;
    static std::shared_ptr<Transposon> ORIGINAL_TE_;

    std::map<uint_fast32_t, std::shared_ptr<Transposon>> sites_;
};

} // namespace tek

#endif /* TEK_HAPLOID_HPP_ */
