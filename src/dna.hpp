// -*- mode: c++; coding: utf-8 -*-
/*! @file dna.hpp
    @brief Interface of DNA class
*/
#pragma once
#ifndef TEK_DNA_HPP_
#define TEK_DNA_HPP_

#include <sfmt.hpp>

#include <cstdint>
#include <array>
#include <iostream>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

template <size_t N>
class DNA {
  public:
    size_t count() const {
        size_t cnt = 0;
        for (const auto x: sequence_) {
            if (x > 0U) ++cnt;
        }
        return cnt;
    }

    template <class URNG> inline
    void flip(const size_t i, URNG& generator) {
        // NOTE: mutation rate is reduced by 1/4
        sequence_[i] ^= (generator() & 0b00000011);
    }

    int_fast8_t operator[](const size_t i) const {
        return sequence_[i];
    }

    std::ostream& write(std::ostream& ost) const {
        for (const auto x: sequence_) {
            ost << static_cast<int_fast16_t>(x);
        }
        return ost;
    }

    friend std::ostream& operator<<(std::ostream& ost, const DNA& x) {
        return x.write(ost);
    }

    static void test() {
        DNA<N> dna;
        std::cerr << dna << std::endl;
        dna.flip(0, wtl::sfmt());
        dna.flip(2, wtl::sfmt());
        dna.flip(4, wtl::sfmt());
        std::cerr << dna << std::endl;
    }

  private:
    std::array<int_fast8_t, N> sequence_ = {};
};

} // namespace tek

#endif /* TEK_DNA_HPP_ */
