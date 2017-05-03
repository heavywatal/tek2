// -*- mode: c++; coding: utf-8 -*-
/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"

#include <iostream>

#include <wtl/debug.hpp>

namespace tek {

Population::Population(const size_t n) {HERE;
    individuals_.reserve(n);
    for (size_t i=0; i<n; ++i) {
        individuals_.push_back(Individual());
    }
}

std::ostream& Population::write(std::ostream& ost) const {HERE;
    for (const auto& x: individuals_) {
        ost << x << std::endl;
    }
    return ost;
}

std::ostream& operator<<(std::ostream& ost, const Population& pop) {
    return pop.write(ost);
}

void Population::unit_test() {HERE;
    Population pop(4);
    std::cout << pop << std::endl;
}

} // namespace tek
