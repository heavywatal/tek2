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
#include <iostream>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

/*! @brief DNA class
*/
template <size_t N>
class DNA {
  public:
    //! default constructor
    DNA() noexcept = default;
    //! default copy constructor
    DNA(const DNA&) = default;
    //! default move constructor
    DNA(DNA&&) noexcept = default;
    //! construct from sequence
    DNA(std::valarray<uint_fast8_t>&& s) noexcept {
        for (uint_fast32_t i=0; i<N; ++i) {
            has_3bonds_.set(i, 0b10u & s[i]);
            is_pyrimidine_.set(i, 0b01u & s[i]);
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
        while ((0b11u & random_bits) == 0u) {
            random_bits >>= 2u;
        }
        if (0b10u & random_bits) {has_3bonds_.flip(i);}
        if (0b01u & random_bits) {is_pyrimidine_.flip(i);}
    }

    //! get i-th nucleotide
    const char& operator[](uint_fast32_t i) const noexcept {
        return translate(get(i));
    }

    //! get i-th nucleotide as integer
    uint_fast8_t get(uint_fast32_t i) const noexcept {
        return (has_3bonds_[i] << 1u) + is_pyrimidine_[i];
    }

    //! translate integer to nucleotide character and print
    std::ostream& write(std::ostream& ost) const {
        for (uint_fast32_t i=0; i<N; ++i) {
            ost << operator[](i);
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
        static constexpr char NUCLEOTIDE[] = "ATGC";
        return NUCLEOTIDE[x];
    }

    // sequence as two bitsets: {00: A, 01: T, 10: G, 11: C}
    //! {0: AT, 1: GC}
    std::bitset<N> has_3bonds_;
    //! {0: AG, 1: TC}
    std::bitset<N> is_pyrimidine_;
};

//! @cond
template <size_t N>
class Homolog {
  public:
    Homolog() noexcept: counts_(N, {0u, 0u, 0u, 0u}) {}

    void collect(const DNA<N>& seq) noexcept {
        for (uint_fast32_t i=0; i<N; ++i) {
            ++counts_[i][seq.get(i)];
        }
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
//! @endcond

} // namespace tek

#endif /* TEK_DNA_HPP_ */
