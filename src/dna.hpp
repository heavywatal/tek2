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

/*! @brief DNA class
*/
template <uint_fast32_t N>
class DNA {
    //! template for translation from integer to character
    static const std::string NUCLEOTIDE;

  public:
    //! diviation from the original
    uint_fast32_t count() const {
        uint_fast32_t cnt = 0;
        for (const auto x: sequence_) {
            if (x > 0U) ++cnt;
        }
        return cnt;
    }

    //! mutate i-th site
    template <class URNG> inline
    void flip(const uint_fast32_t i, URNG& generator) {
        typename URNG::result_type random_bits = 0U;
        while ((random_bits = generator()) < 1U) {;}
        uint_fast8_t two_bits = 0U;
        while ((two_bits = random_bits & 0b00000011) < 1U) {
            random_bits >>= 2;
        }
        sequence_[i] ^= two_bits;
    }

    //! get i-th nucleotide as char
    const char& operator[](const uint_fast32_t i) const {
        return NUCLEOTIDE[sequence_[i]];
    }

    //! translate integer to nucleotide character and print
    std::ostream& write(std::ostream& ost) const {
        for (const auto x: sequence_) {
            ost << NUCLEOTIDE[x];
        }
        return ost;
    }

    //! shortcut of write()
    friend std::ostream& operator<<(std::ostream& ost, const DNA& x) {
        return x.write(ost);
    }

    //! count different nucleotides
    uint_fast32_t operator-(const DNA<N>& other) const {
        uint_fast32_t diff = 0;
        for (uint_fast32_t i=0; i<N; ++i) {
            if (this->sequence_[i] != other.sequence_[i]) ++diff;
        }
        return diff;
    }

  private:
    //! sequence as integer array
    std::array<uint_fast8_t, N> sequence_ = {};
};

template <uint_fast32_t N>
const std::string DNA<N>::NUCLEOTIDE = "ACGT";

//! unit test
template <class URNG> inline
void DNA_test(URNG& generator) {
    constexpr uint_fast32_t N = 30;
    DNA<N> x, y;
    std::cerr << x << std::endl;
    for (uint_fast32_t i=0; i<N; ++i) {
        x.flip(i, generator);
    }
    std::cerr << x << std::endl;
    std::cerr << x - y << std::endl;
}

} // namespace tek

#endif /* TEK_DNA_HPP_ */
