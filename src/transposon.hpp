// -*- mode: c++; coding: utf-8 -*-
/*! @file transposon.hpp
    @brief Interface of Transposon class
*/
#pragma once
#ifndef TEK_TRANSPOSON_HPP_
#define TEK_TRANSPOSON_HPP_

#include <iosfwd>
#include <bitset>

#include <boost/program_options.hpp>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

class Transposon {
  public:
    static constexpr size_t LENGTH = 300;

    Transposon() = default;

    void mutate();
    void indel() {has_indel_ = true;}
    double activity() const;
    double transposition_rate() const {
        return MAX_TRANSPOSITION_RATE_ * activity();
    }

    double dn() const {return nonsynonymous_sites_.count() * OVER_NONSYNONYMOUS_SITES_;}
    double ds() const {return synonymous_sites_.count() * OVER_SYNONYMOUS_SITES_;}

    std::ostream& write_summary(std::ostream&) const;
    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Transposon&);

    static void set_parameters();
    static boost::program_options::options_description options_desc();
    static void unit_test();

  private:

    static constexpr double MAX_TRANSPOSITION_RATE_ = 0.01;
    static constexpr size_t NUM_SYNONYMOUS_SITES_ = LENGTH / 3;
    static constexpr size_t NUM_NONSYNONYMOUS_SITES_ = LENGTH - NUM_SYNONYMOUS_SITES_;
    static constexpr double OVER_SYNONYMOUS_SITES_ = 1.0 / NUM_SYNONYMOUS_SITES_;
    static constexpr double OVER_NONSYNONYMOUS_SITES_ = 1.0 / NUM_NONSYNONYMOUS_SITES_;
    static double ALPHA_;
    static double THRESHOLD_;
    static unsigned int BETA_;

    std::bitset<NUM_NONSYNONYMOUS_SITES_> nonsynonymous_sites_;
    std::bitset<LENGTH - NUM_NONSYNONYMOUS_SITES_> synonymous_sites_;
    bool has_indel_ = false;
};

} // namespace tek

#endif /* TEK_TRANSPOSON_HPP_ */
