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
class DNA {
  public:
    //! construct
    DNA(size_t n): sequence_(n) {}
    //! construct from sequence
    DNA(std::valarray<uint_fast8_t>&& s)
    : sequence_(std::forward<std::valarray<uint_fast8_t>>(s)) {}

    //! diviation from the original
    uint_fast32_t count() const {
        return static_cast<uint_fast32_t>(wtl::count(sequence_ > 0u));
    }

    //! mutate i-th site
    template <class URBG> inline
    void flip(size_t i, URBG& engine) {
        typename URBG::result_type random_bits = 0u;
        while ((random_bits = engine()) == 0u) {;}
        while ((0b00000011u & random_bits) == 0u) {
            random_bits >>= 2;
        }
        sequence_[i] ^= (0b00000011u & random_bits);
    }

    //! get i-th nucleotide
    uint_fast8_t operator[](size_t i) const {
        return sequence_[i];
    }

    //! get i-th nucleotide as char
    const char& at(size_t i) const {
        return translate(sequence_[i]);
    }

    //! translate integer to nucleotide character and print
    std::ostream& write(std::ostream& ost) const {
        for (auto x: sequence_) {
            ost << translate(x);
        }
        return ost;
    }

    //! shortcut of write()
    friend std::ostream& operator<<(std::ostream& ost, const DNA& x) {
        return x.write(ost);
    }

    //! comparison of sequence
    std::valarray<bool> operator!=(const DNA& other) const {
        return this->sequence_ != other.sequence_;
    }

  private:
    //! translate integer to character
    static const char& translate(uint_fast8_t x) {
        static const std::string NUCLEOTIDE = "ACGT";
        return NUCLEOTIDE[x];
    }

    //! sequence as integer array
    std::valarray<uint_fast8_t> sequence_;
};

/*! @brief Homolog class
*/
class Homolog {
  public:
    Homolog(size_t length): counts_(length, {0u, 0u, 0u, 0u}) {}

    Homolog& operator+=(const DNA& seq) {
        const size_t n = counts_.size();
        for (size_t i=0; i<n; ++i) {
            ++counts_[i][seq[i]];
        }
        return *this;
    }

    DNA majority() const {
        const size_t n = counts_.size();
        std::valarray<uint_fast8_t> result(n);
        for (size_t i=0; i<n; ++i) {
            const auto& v = counts_[i];
            result[i] = static_cast<uint_fast8_t>(std::distance(v.begin(), std::max_element(v.begin(), v.end())));
        }
        return DNA(std::move(result));
    }

  private:
    std::vector<std::vector<uint_fast32_t>> counts_;
};

} // namespace tek

#endif /* TEK_DNA_HPP_ */
