// -*- mode: c++; coding: utf-8 -*-
/*! @file program.hpp
    @brief Interface of Program class
*/
#pragma once
#ifndef TEK_PROGRAM_HPP_
#define TEK_PROGRAM_HPP_

#include <vector>
#include <string>

#include <boost/program_options.hpp>

namespace tek {

/*! @brief Represents single run
*/
class Program {
  public:
    //! Parse command arguments
    Program(const std::vector<std::string>& args);
    //! Parse command arguments
    Program(int argc, char* argv[])
    : Program(std::vector<std::string>(argv, argv + argc)) {}

    //! Top level function that should be called once from global main
    void run();

  private:
    //! called from run()
    void main();
    //! options description for Program class
    boost::program_options::options_description options_desc();
    //! options description for positional arguments
    boost::program_options::options_description positional_desc();
    //! Print help message and exit
    void help_and_exit();

    //! population size
    size_t popsize_ = 500;
    //! initial number of individuals with TE
    size_t initial_freq_ = 1;
    //! maximum number of generations to simulate
    size_t num_generations_ = 1000;
    //! number of generations to simulate after population split
    size_t num_generations_after_split_ = 0;
    //! interval of recording
    size_t record_interval_ = 10;
    //! enum Recording
    int record_flags_ = 3;
    //! number of threads
    unsigned int concurrency_ = 1;
    //! dummy
    size_t num_repeats_ = 1;
    //! name of output directory
    std::string outdir_ = "";
    //! writen to "program_options.conf"
    std::string config_string_;
};

} // namespace tek

#endif /* TEK_PROGRAM_HPP_ */
