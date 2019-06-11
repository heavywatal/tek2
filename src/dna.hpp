/*! @file dna.hpp
    @brief Interface of DNA class
*/
#pragma once
#ifndef TEK_DNA_HPP_
#define TEK_DNA_HPP_

#include <cstdint>
#include <bitset>
#include <valarray>
#include <vector>
#include <numeric>
#include <string>
#include <iostream>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

/*! @brief DNA class
*/
template <size_t N>
class DNA {
  public:
    //! default constructor
    DNA() = default;
    //! default copy constructor
    DNA(const DNA&) = default;
    //! move constructor
    DNA(DNA&&) noexcept = default;
    // //! construct from sequence
    DNA(std::valarray<uint_fast8_t>&& s) noexcept {
        for (unsigned i=0; i<N; ++i) {
            is_pyrimidine_.set(i, (0b00000010u & s[i]) >> 1);
            has_3bonds_.set(i, 0b00000001u & s[i]);
        }
    }

    //! diviation from the original
    uint_fast32_t count() const noexcept {
        return (is_pyrimidine_ | has_3bonds_).count();
    }

    //! mutate i-th site
    template <class URBG> inline
    void flip(uint_fast32_t i, URBG& engine) noexcept {
        typename URBG::result_type random_bits = 0u;
        while ((random_bits = engine()) == 0u) {;}
        while ((0b00000011u & random_bits) == 0u) {
            random_bits >>= 2;
        }
        is_pyrimidine_.set(i, (0b00000010u & random_bits) >> 1);
        has_3bonds_.set(i, 0b00000001u & random_bits);
    }

    //! get i-th nucleotide
    uint_fast8_t operator[](uint_fast32_t i) const noexcept {
        return (is_pyrimidine_[i] << 1) + has_3bonds_[i];
    }
    //! get i-th nucleotide as char
    const char& at(uint_fast32_t i) const noexcept {
        return translate(operator[](i));
    }

    //! translate integer to nucleotide character and print
    std::ostream& write(std::ostream& ost) const {
        for (unsigned int i=0; i<N; ++i) {
            ost << at(i);
        }
        return ost;
    }

    //! shortcut of write()
    friend std::ostream& operator<<(std::ostream& ost, const DNA& x) {
        return x.write(ost);
    }

    //! Hamming distance
    uint_fast32_t operator-(const DNA& other) const noexcept {
        return ((this->is_pyrimidine_ ^ other.is_pyrimidine_) | (this->has_3bonds_ ^ other.has_3bonds_)).count();
    }

  private:
    //! translate integer to character
    static const char& translate(uint_fast8_t x) noexcept {
        static const std::string NUCLEOTIDE = "ACGT";
        return NUCLEOTIDE[x];
    }
    //! count the number of trues
    static uint_fast32_t count(const std::valarray<bool>& v) noexcept {
        return std::accumulate(std::begin(v), std::end(v), 0u,
          [](uint_fast32_t x, bool b) {
              if (b) {return ++x;} else {return x;}
          });
    }

    //! sequence as two bitsets: {00: A, 01: G, 10: T, 11: C}
    std::bitset<N> is_pyrimidine_; // {0: AG, 1: TC}
    std::bitset<N> has_3bonds_;    // {0: AT, 1: GC}
};

/*! @brief Homolog class
*/
template <size_t N>
class Homolog {
  public:
    Homolog() noexcept: counts_(N, {0u, 0u, 0u, 0u}) {}

    Homolog& operator+=(const DNA<N>& seq) noexcept {
        for (uint_fast32_t i=0; i<N; ++i) {
            ++counts_[i][seq[i]];
        }
        return *this;
    }

    DNA<N> majority() const noexcept {
        std::valarray<uint_fast8_t> result(N);
        for (uint_fast32_t i=0; i<N; ++i) {
            const auto& v = counts_[i];
            result[i] = static_cast<uint_fast8_t>(std::distance(v.begin(), std::max_element(v.begin(), v.end())));
        }
        return DNA<N>(std::move(result));
    }

  private:
    std::vector<std::vector<uint_fast32_t>> counts_;
};

} // namespace tek

#endif /* TEK_DNA_HPP_ */
