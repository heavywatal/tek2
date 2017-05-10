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
    Transposon() = default;

    void mutate();
    void indel() {has_indel_ = true;}
    double transposition_rate() const {
        return MAX_TRANSPOSITION_RATE_ * activity();
    }

    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Transposon&);

    static void set_parameters();
    static boost::program_options::options_description options_desc();
    static void unit_test();

  private:
    double activity() const;

    static constexpr double MAX_TRANSPOSITION_RATE_ = 0.01;
    static constexpr size_t LENGTH_ = 300;
    static constexpr size_t NUM_NONSYNONYMOUS_SITES_ = LENGTH_ * 2 / 3;
    static constexpr double OVER_NONSYNONYMOUS_SITES_ = 1.0 / NUM_NONSYNONYMOUS_SITES_;
    static double ALPHA_;
    static double THRESHOLD_;
    static unsigned int BETA_;

    std::bitset<NUM_NONSYNONYMOUS_SITES_> nonsynonymous_sites_;
    std::bitset<LENGTH_ - NUM_NONSYNONYMOUS_SITES_> synonymous_sites_;
    bool has_indel_ = false;
};

} // namespace tek

#endif /* TEK_TRANSPOSON_HPP_ */
