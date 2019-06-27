/*! @file population.hpp
    @brief Interface of Population class
*/
#pragma once
#ifndef TEK_POPULATION_HPP_
#define TEK_POPULATION_HPP_

#include <iosfwd>
#include <vector>
#include <random>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tek {

class Haploid;

//! bits to denote what to record
enum class Recording: int {
    none     = 0b00000000,
    activity = 0b00000001,
    sequence = 0b00000010,
    fitness  = 0b00000100,
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

//! @brief Parameters for Population class
/*! @ingroup params
*/
struct PopulationParams {
    //! number of individuals to sample
    size_t SAMPLE_SIZE = 10u;
    //! number of threads
    unsigned int CONCURRENCY = 1u;
    //! max number of species that can coexist at a time
    unsigned int MAX_COEXISTENCE = 42u;
};

/*! @brief Population class
*/
class Population {
  public:
    //! Alias
    using param_type = PopulationParams;

    //! @addtogroup params
    //! @{

    //! \f$\theta = 4N\mu\f$, population mutation rate
    static constexpr double THETA = 0.01;
    //! \f$\rho = 4Nc\f$, population recombination rate
    static constexpr double RHO = 20000;
    //! @} params /2/////////3/////////4/////////5/////////6/////////7/////////

    //! constructor
    Population(size_t size, size_t num_founders=1);
    //! default copy constructor
    Population(const Population& other) = default;
    //! destructor
    ~Population();

    //! return false if TE is extinct
    bool evolve(size_t max_generations, size_t record_interval,
                Recording flags=Recording::activity | Recording::fitness,
                size_t t_hyperactivate = 0u);

    //! write summary in JSON format
    std::ostream& write_summary(std::ostream&) const;
    //! count identicals and write FASTA for i-th individual
    std::ostream& write_fasta_individual(std::ostream&, size_t i) const;
    //! call write_fasta_individual() repeatedly
    std::ostream& write_fasta(std::ostream&, size_t num_individuals=-1u) const;
    friend std::ostream& operator<<(std::ostream&, const Population&);

    //! Set #PARAM_
    static void param(const param_type& p) {PARAM_ = p;}
    //! Get #PARAM_
    static const param_type& param() {return PARAM_;}
    //! Set #SEEDER_ seed
    static void seed(std::mt19937_64::result_type value) {SEEDER_.seed(value);}

  private:
    //! Parameters shared among instances
    static param_type PARAM_;
    //! seed generator for Haploid::URBG
    static std::mt19937_64 SEEDER_;

    //! proceed one generation and return fitness record
    std::vector<double> step(double previous_max_fitness=1.0);
    //! find farthest element, count species, and cause speciation if qualified
    void eval_species_distance();
    //! return true if no TE exists in #gametes_
    bool is_extinct() const;
    //! summarize and write activity
    void write_activity(std::ostream&, size_t time, bool header) const;

    //! vector of chromosomes, not individuals
    std::vector<Haploid> gametes_;
};

} // namespace tek

#endif /* TEK_POPULATION_HPP_ */
