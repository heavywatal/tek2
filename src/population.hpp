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

//! bits to denote what to record
enum class Recording: int {
    none     = 0b00000000,
    activity = 0b00000001,
    fitness  = 0b00000010,
    sequence = 0b00000100,
};

//! operator OR
constexpr Recording operator&(Recording x, Recording y) {
    return static_cast<Recording>(static_cast<int>(x) & static_cast<int>(y));
}

//! operator AND
constexpr Recording operator|(Recording x, Recording y) {
    return static_cast<Recording>(static_cast<int>(x) | static_cast<int>(y));
}

class Haploid;

/*! @brief Population class
*/
class Population {
  public:
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    //! @addtogroup params
    //! @{

    //! \f$\theta = 4N\mu\f$, population mutation rate
    static constexpr double THETA = 0.01;
    //! \f$\rho = 4Nc\f$, population recombination rate
    static constexpr double RHO = 200;
    //! @} params
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! constructor
    Population(const size_t size, const size_t num_founders=1, const unsigned int concurrency=1);
    //! default copy constructor
    Population(const Population& other) = default;

    //! return false if TE is extinct
    bool evolve(const size_t max_generations, const size_t record_interval,
                const Recording flags=Recording::activity | Recording::fitness);

    //! write summary in JSON format
    std::ostream& write_summary(std::ostream&) const;
    //! count identicals and write FASTA
    std::ostream& write_fasta(std::ostream&) const;
    //! write_fasta() for i-th individual
    std::ostream& write_individual(std::ostream&, const size_t i=0) const;
    friend std::ostream& operator<<(std::ostream&, const Population&);
    //! unit test
    static void test();
  private:
    //! proceed one generation and return fitness record
    std::vector<double> step(const double previous_max_fitness=1.0);
    //! return true if no TE exists in #gametes_
    bool is_extinct() const;
    //! vector of chromosomes, not individuals
    std::vector<Haploid> gametes_;
    //! number of threads
    const unsigned int concurrency_;
};

} // namespace tek

#endif /* TEK_POPULATION_HPP_ */
