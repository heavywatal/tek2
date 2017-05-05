// -*- mode: c++; coding: utf-8 -*-
/*! @file individual.hpp
    @brief Interface of Individual class
*/
#pragma once
#ifndef TEK_INDIVIDUAL_HPP_
#define TEK_INDIVIDUAL_HPP_

#include <iosfwd>

#include "transposon.hpp"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

class Individual {
  public:
    Individual() = default;

    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Individual&);

    static void unit_test();

  private:
    double fitness_ = 1.0;
};

} // namespace tek

#endif /* TEK_INDIVIDUAL_HPP_ */
