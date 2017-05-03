// -*- mode: c++; coding: utf-8 -*-
/*! @file population.hpp
    @brief Interface of Population class
*/
#pragma once
#ifndef TEK_POPULATION_HPP_
#define TEK_POPULATION_HPP_

#include <vector>

#include "individual.hpp"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

class Population {
  public:
    Population(const size_t n);

    static void unit_test() {}
  private:
    std::vector<Individual> individuals_;
};

} // namespace tek

#endif /* TEK_POPULATION_HPP_ */
