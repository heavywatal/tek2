// -*- mode: c++; coding: utf-8 -*-
/*! @file population.hpp
    @brief Interface of Population class
*/
#pragma once
#ifndef TEK_POPULATION_HPP_
#define TEK_POPULATION_HPP_

#include <iosfwd>
#include <vector>
#include <map>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

class Haploid;

class Population {
  public:
    static constexpr double THETA = 0.01;
    static constexpr double RHO = 200;

    Population(const size_t size, const size_t num_founders=1, const unsigned int concurrency=1);
    bool evolve(const size_t max_generations);
    void sample() const;

    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Population&);

    static void test();
  private:
    void step();

    bool is_extinct() const;

    std::vector<Haploid> gametes_;
    const unsigned int concurrency_;
};

} // namespace tek

#endif /* TEK_POPULATION_HPP_ */
