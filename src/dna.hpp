/*! @file dna.hpp
    @brief Interface of DNA class
*/
#pragma once
#ifndef TEK_DNA_HPP_
#define TEK_DNA_HPP_

#include <wtl/numeric.hpp>

#include <cstdint>
#include <valarray>
#include <string>
#include <iostream>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

/*! @brief DNA class
*/
template <size_t N>
class DNA {
    //! template for translation from integer to character
    static const std::string NUCLEOTIDE;

  public:
    DNA(): sequence_(N) {}

    //! diviation from the original
    uint_fast32_t count() const {
        return wtl::count(sequence_ > 0u);
    }

    //! mutate i-th site
    template <class URBG> inline
    void flip(const size_t i, URBG& generator) {
        typename URBG::result_type random_bits = 0u;
        while ((random_bits = generator()) == 0u) {;}
        while ((0b00000011u & random_bits) == 0u) {
            random_bits >>= 2;
        }
        sequence_[i] ^= (0b00000011u & random_bits);
    }

    //! get i-th nucleotide as char
    const char& operator[](const size_t i) const {
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

    //! comparison of sequence
    std::valarray<bool> operator!=(const DNA<N>& other) const {
        return this->sequence_ != other.sequence_;
    }

  private:
    //! sequence as integer array
    std::valarray<uint_fast8_t> sequence_;
};

template <size_t N>
const std::string DNA<N>::NUCLEOTIDE = "ACGT";

//! unit test
template <class URBG> inline
void DNA_test(URBG& generator) {
    constexpr size_t N = 30u;
    DNA<N> x, y;
    std::cerr << x << std::endl;
    for (size_t i=0u; i<N; ++i) {
        x.flip(i, generator);
    }
    std::cerr << x << std::endl;
    std::cerr << wtl::count(x != y) << std::endl;
}

} // namespace tek

#endif /* TEK_DNA_HPP_ */
