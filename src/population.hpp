/*! @file population.hpp
    @brief Interface of Population class
*/
#pragma once
#ifndef TEK_POPULATION_HPP_
#define TEK_POPULATION_HPP_

#include "haploid.hpp"

#include <iosfwd>
#include <vector>

namespace boost {namespace program_options {class options_description;}}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

//! bits to denote what to record
enum class Recording: int {
    none     = 0b00000000,
    activity = 0b00000001,
    fitness  = 0b00000010,
    sequence = 0b00000100,
    summary  = 0b00001000,
};

//! operator OR
constexpr Recording operator&(Recording x, Recording y) {
    return static_cast<Recording>(static_cast<int>(x) & static_cast<int>(y));
}

//! operator AND
constexpr Recording operator|(Recording x, Recording y) {
    return static_cast<Recording>(static_cast<int>(x) | static_cast<int>(y));
}

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
    static constexpr double RHO = 20000;
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
    //! count identicals and write FASTA for i-th individual
    std::ostream& write_fasta_individual(std::ostream&, size_t i) const;
    //! call write_fasta_individual() repeatedly
    std::ostream& write_fasta(std::ostream&, size_t num_individuals=-1u) const;
    friend std::ostream& operator<<(std::ostream&, const Population&);

    //! options description for Population class
    static boost::program_options::options_description options_desc();

  private:
    //! proceed one generation and return fitness record
    std::vector<double> step(const double previous_max_fitness=1.0);
    //! find farthest element, count species, and cause speciation if qualified
    void supply_new_species();
    //! return true if no TE exists in #gametes_
    bool is_extinct() const;
    //! summarize and write activity
    void write_activity(std::ostream&, const size_t time, const bool header) const;

    //! number of individuals to sample
    static size_t SAMPLE_SIZE_;

    //! vector of chromosomes, not individuals
    std::vector<Haploid> gametes_;
    //! number of threads
    const unsigned int concurrency_;
};

} // namespace tek

#endif /* TEK_POPULATION_HPP_ */
