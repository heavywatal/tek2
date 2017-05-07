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
    static constexpr unsigned int SEQUENCE_LENGTH_ = 200;
    static constexpr double OVER_L_ = 1.0 / SEQUENCE_LENGTH_;
    static double ALPHA_;
    static double THRESHOLD_;
    static unsigned int BETA_;

    std::bitset<SEQUENCE_LENGTH_> sequence_;
    bool has_indel_ = false;
};

} // namespace tek

#endif /* TEK_TRANSPOSON_HPP_ */
