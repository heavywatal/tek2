/*! @file program.hpp
    @brief Interface of Program class
*/
#pragma once
#ifndef TEK_PROGRAM_HPP_
#define TEK_PROGRAM_HPP_

#include <vector>
#include <string>

namespace tek {

/*! @brief Program class
*/
class Program {
  public:
    //! Parse command arguments
    Program(const std::vector<std::string>& args);
    //! Parse command arguments
    Program(int argc, char* argv[])
    : Program(std::vector<std::string>(argv + 1, argv + argc)) {}

    //! Top level function that should be called once from global main
    void run();

  private:
    //! called from run()
    void main();

    //! writen to "config.json"
    std::string config_;
};

} // namespace tek

#endif /* TEK_PROGRAM_HPP_ */
