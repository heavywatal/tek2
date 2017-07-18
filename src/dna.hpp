// -*- mode: c++; coding: utf-8 -*-
/*! @file dna.hpp
    @brief Interface of DNA class
*/
#pragma once
#ifndef TEK_DNA_HPP_
#define TEK_DNA_HPP_

#include <cstdint>
#include <array>
#include <string>
#include <iostream>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

template <size_t N>
class DNA {
    static const std::string NUCLEOTIDE;

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

    uint_fast8_t operator[](const size_t i) const {
        return sequence_[i];
    }

    char& nucleotide(const size_t i) const {
        return NUCLEOTIDE[sequence_[i]];
    }

    std::ostream& write(std::ostream& ost) const {
        for (const auto x: sequence_) {
            ost << NUCLEOTIDE[x];
        }
        return ost;
    }

    friend std::ostream& operator<<(std::ostream& ost, const DNA& x) {
        return x.write(ost);
    }

    template <class URNG> inline
    static void test(URNG& generator) {
        DNA<N> dna;
        std::cerr << dna << std::endl;
        for (size_t i=0; i<N; ++i) {
            dna.flip(i, generator);
        }
        std::cerr << dna << std::endl;
    }

  private:
    std::array<uint_fast8_t, N> sequence_ = {};
};

template <size_t N>
const std::string DNA<N>::NUCLEOTIDE = "ACGT";

} // namespace tek

#endif /* TEK_DNA_HPP_ */
