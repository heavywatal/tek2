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

class Haploid {
  public:
    using URNG = wtl::sfmt19937;

    Haploid() = default;
    Haploid(Haploid&& other) = default;
    Haploid& operator=(Haploid&&) = default;

    Haploid gametogenesis(const Haploid& other, URNG& rng) const;
    void transpose_mutate(Haploid& other, URNG& rng);
    double fitness(const Haploid&) const;

    bool has_transposon() const {return copy_number_ > 0U;};
    std::map<double, unsigned int> count_activity() const;

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
    Haploid(const Haploid&) = default;
    Haploid& operator=(const Haploid&) = default;

    static void set_SELECTION_COEFS_GP();

    std::vector<std::shared_ptr<Transposon>> transpose(URNG&);
    void recombine(Haploid&, URNG&);
    void mutate(URNG&);
    void evaluate_sites();

    static void test_selection_coefs_cn();
    static void test_selection_coefs_gp();
    static void test_recombination();

    static constexpr uint_fast32_t NUM_SITES = 2'000'000;
    static constexpr double INDEL_RATIO_ = 0.2;
    static constexpr double TAU_ = 1.5;
    static constexpr double PROP_FUNCTIONAL_SITES_ = 0.75;
    static double XI_;
    static double EXCISION_RATE_;
    static double MEAN_SELECTION_COEF_;

    static double RECOMBINATION_RATE_;
    static double INDEL_RATE_;
    static std::valarray<double> SELECTION_COEFS_GP_;
    static std::uniform_int_distribution<uint_fast32_t> UNIFORM_SITES_;
    static std::poisson_distribution<unsigned int> NUM_MUTATIONS_DIST_;
    static std::shared_ptr<Transposon> ORIGINAL_TE_;

    std::map<uint_fast32_t, std::shared_ptr<Transposon>> sites_;
    double prod_1_zs_ = 1.0;
    uint_fast32_t copy_number_ = 0;
};

} // namespace tek

#endif /* TEK_HAPLOID_HPP_ */
