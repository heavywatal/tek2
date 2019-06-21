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
#include <sstream>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

template <size_t N> inline
std::bitset<2 * N> make_odd_bits() {
    std::ostringstream oss;
    for (size_t i=0; i<N; ++i) {
        oss << "01";
    }
    return std::bitset<2 * N>(oss.str());
}

template <size_t N> inline
std::bitset<N> odd_and(const std::bitset<N>& x) {
    static const std::bitset<N> mask = make_odd_bits<N / 2>();
    return mask & x;
}

/*! @brief DNA class
*/
template <size_t N>
class DNA {
  public:
    //! default constructor
    DNA() = default;
    //! default copy constructor
    DNA(const DNA&) = default;
    //! default move constructor
    DNA(DNA&&) noexcept = default;
    //! construct from sequence
    DNA(std::valarray<uint_fast8_t>&& s) noexcept {
        for (unsigned i=0; i<N; ++i) {
            bits_ |= s[N - i];
            bits_ <<= 2u;
            // bits_.set(i + N, 0b10u & s[i]);
            // bits_.set(i, 0b01u & s[i]);
            // has_3bonds_.set(i, 0b10u & s[i]);
            // is_pyrimidine_.set(i, 0b01u & s[i]);
        }
    }

    //! mutate i-th site
    template <class URBG> inline
    void flip(uint_fast32_t i, URBG& engine) noexcept {
        typename URBG::result_type random_bits = 0u;
        while ((random_bits = engine()) == 0u) {;}
        while ((0b11u & random_bits) == 0u) {
            random_bits >>= 2u;
        }
        std::bitset<2 * N> mask(0b00000011u & random_bits);
        mask <<= (2u * i);
        bits_ ^= mask;
        // if (0b10u & random_bits) {bits_.flip(i + N);}
        // if (0b01u & random_bits) {bits_.flip(i)};
        // if (0b10u & random_bits) {has_3bonds_.flip(i);}
        // if (0b01u & random_bits) {is_pyrimidine_.flip(i);}
    }

    //! get i-th nucleotide
    const char& operator[](uint_fast32_t i) const noexcept {
        return translate(get(i));
    }

    //! get i-th nucleotide as integer
    uint_fast8_t get(uint_fast32_t i) const noexcept {
        static const std::bitset<2 * N> mask(0b00000011u);
        return ((bits_ >> (2u * i)) & mask).to_ulong();
        // return (bits_[i + N] << 1) + bits_[i];
        // return (has_3bonds_[i] << 1) + is_pyrimidine_[i];
    }

    //! translate integer to nucleotide character and print
    std::ostream& write(std::ostream& ost) const {
        for (unsigned int i=0; i<N; ++i) {
            ost << operator[](i);
        }
        return ost;
    }

    //! shortcut of write()
    friend std::ostream& operator<<(std::ostream& ost, const DNA& x) {
        return x.write(ost);
    }

    //! diviation from the original
    uint_fast32_t count() const noexcept {
        return (odd_and(bits_) | odd_and(bits_ >> 1u)).count();
        // return ((bits_ >> N) | (bits_ << N >> N)).count();
        // return (is_pyrimidine_ | has_3bonds_).count();
    }

    //! Hamming distance
    uint_fast32_t operator-(const DNA& other) const noexcept {
        auto diff = bits_ ^ other.bits_;
        return (odd_and(diff) | odd_and(diff >> 1u)).count();
        // return ((diff >> N) | (diff << N >> N)).count();
        // return ((this->is_pyrimidine_ ^ other.is_pyrimidine_) | (this->has_3bonds_ ^ other.has_3bonds_)).count();
    }

  private:
    //! translate integer to character
    static const char& translate(uint_fast8_t x) noexcept {
        static const std::string NUCLEOTIDE = "ATGC";
        return NUCLEOTIDE[x];
    }
    //! count the number of trues
    static uint_fast32_t count(const std::valarray<bool>& v) noexcept {
        return std::accumulate(std::begin(v), std::end(v), 0u,
          [](uint_fast32_t x, bool b) {
              if (b) {return ++x;} else {return x;}
          });
    }

    //! sequence as two bitsets: {00: A, 01: T, 10: G, 11: C}
    std::bitset<2 * N> bits_;
};

/*! @brief Homolog class
*/
template <size_t N>
class Homolog {
  public:
    Homolog() noexcept: counts_(N, {0u, 0u, 0u, 0u}) {}

    Homolog& operator+=(const DNA<N>& seq) noexcept {
        for (uint_fast32_t i=0; i<N; ++i) {
            ++counts_[i][seq.get(i)];
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
