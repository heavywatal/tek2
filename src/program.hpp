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
    Program(int argc, char* argv[])
    : Program(std::vector<std::string>(argv, argv + argc)) {}

    //! Top level function that should be called once from main()
    void run();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    boost::program_options::options_description options_desc();
    boost::program_options::options_description positional_desc();
    void help_and_exit();

    size_t popsize_ = 500;
    size_t initial_freq_ = 1;
    size_t num_generations_ = 1000;
    size_t record_interval_ = 10;
    int record_flags_ = 3;
    unsigned int concurrency_ = 1;
    size_t num_repeats_ = 1;
    std::string outdir_ = "";
};

} // namespace tek

#endif /* TEK_PROGRAM_HPP_ */
