// -*- mode: c++; coding: utf-8 -*-
/*! @file population.cpp
    @brief Implementation of Population class
*/

#include "population.hpp"

namespace tek {

Population::Population(const size_t n) {
    individuals_.reserve(n);
    for (size_t i=0; i<n; ++i) {
        individuals_.push_back(Individual());
    }
}

} // namespace tek
