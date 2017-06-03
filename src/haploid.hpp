// -*- mode: c++; coding: utf-8 -*-
/*! @file haploid.hpp
    @brief Interface of Haploid class
*/
#pragma once
#ifndef TEK_HAPLOID_HPP_
#define TEK_HAPLOID_HPP_

#include <boost/program_options.hpp>

#include <iosfwd>
#include <vector>
#include <valarray>
#include <map>
#include <memory>
#include <random>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

class Transposon;

class Haploid {
  public:
    static constexpr unsigned int NUM_SITES = 2000;

    Haploid(): sites_(NUM_SITES) {}
    Haploid(Haploid&& other) = default;
    Haploid& operator=(Haploid&&) = default;

    std::pair<Haploid, Haploid> gametogenesis(const Haploid& other) const {
        Haploid lhalf(*this), rhalf(other);
        lhalf.transpose(rhalf);
        lhalf.recombine(rhalf);
        lhalf.mutate();
        rhalf.mutate();
        return std::make_pair(std::move(lhalf), std::move(rhalf));
    }

    void count_activities(std::map<double, unsigned int>* const) const;

    double fitness(const Haploid&) const;
    unsigned int count_transposons() const;
    bool has_transposon() const;

    std::ostream& write_sample(std::ostream&) const;
    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Haploid&);

    static Haploid copy_founder();
    static void set_parameters(const size_t popsize, const double theta, const double rho);
    static boost::program_options::options_description options_desc();
    static void unit_test();

  private:
    Haploid(const Haploid&) = default;
    Haploid& operator=(const Haploid&) = default;

    static size_t random_index();
    static void set_SELECTION_COEFS_GP();

    std::vector<std::shared_ptr<Transposon>> transpose();
    void transpose(Haploid&);
    void recombine(Haploid&);
    void mutate();

    static constexpr double INDEL_RATIO_ = 0.2;
    static constexpr double TAU_ = 1.5;
    static constexpr double PROP_FUNCTIONAL_SITES_ = 0.75;
    static double XI_;
    static double EXCISION_RATE_;
    static double MEAN_SELECTION_COEF_;
    static std::valarray<double> SELECTION_COEFS_GP_;
    static std::bernoulli_distribution EXCISION_DIST_;
    static std::bernoulli_distribution CHIASMA_DIST_;
    static std::poisson_distribution<> NUM_MUTATIONS_DIST_;
    static std::bernoulli_distribution INDEL_DIST_;
    static std::shared_ptr<Transposon> ORIGINAL_TE_;

    std::vector<std::shared_ptr<Transposon>> sites_;
};

} // namespace tek

#endif /* TEK_HAPLOID_HPP_ */
