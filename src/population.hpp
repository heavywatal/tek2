// -*- mode: c++; coding: utf-8 -*-
/*! @file population.hpp
    @brief Interface of Population class
*/
#pragma once
#ifndef TEK_POPULATION_HPP_
#define TEK_POPULATION_HPP_

#include <vector>

#include "haploid.hpp"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

class Population {
  public:
    Population(const size_t n);

    void step();

    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Population&);

    static void unit_test();
  private:
    std::vector<Haploid> gametes_;
};

} // namespace tek

#endif /* TEK_POPULATION_HPP_ */
