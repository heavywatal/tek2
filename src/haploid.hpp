// -*- mode: c++; coding: utf-8 -*-
/*! @file haploid.hpp
    @brief Interface of Haploid class
*/
#pragma once
#ifndef TEK_HAPLOID_HPP_
#define TEK_HAPLOID_HPP_

#include <iosfwd>
#include <vector>
#include <valarray>
#include <memory>

#include "transposon.hpp"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

class Haploid {
  public:
    static constexpr unsigned int NUM_SITES = 2000;

    Haploid();

    double fitness(const Haploid&) const;
    void mutate(Haploid& other) {
        transpose(other);
        recombine(other);
        mutate();
        other.mutate();
    }

    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Haploid&);

    static void unit_test();

  private:
    std::valarray<double> valarray() const {
        std::valarray<double> v(NUM_SITES);
        for (size_t j=0; j<NUM_SITES; ++j) {
            if (sites_[j]) ++v[j];
        }
        return v;
    }
    static size_t random_index();
    static double s_cn(const unsigned int);

    void transpose(Haploid&);
    std::vector<std::shared_ptr<Transposon>> transpose();
    void recombine(Haploid&);
    void mutate();

    static constexpr double INDEL_RATIO_ = 0.2;
    static constexpr double TAU_ = 1.5;
    static double XI_;
    static double EXCISION_RATE_;
    static std::uniform_int_distribution<size_t> SITES_DIST_;
    static std::shared_ptr<Transposon> ORIGINAL_TE_;

    std::vector<std::shared_ptr<Transposon>> sites_;
};

} // namespace tek

#endif /* TEK_HAPLOID_HPP_ */
