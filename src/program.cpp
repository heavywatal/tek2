// -*- mode: c++; coding: utf-8 -*-
/*! @file program.cpp
    @brief Implementation of Program class
*/
#include "program.hpp"
#include "population.hpp"
#include "haploid.hpp"
#include "transposon.hpp"

#include <wtl/exception.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/os.hpp>
#include <wtl/getopt.hpp>

#include <iostream>

namespace tek {

namespace po = boost::program_options;

inline po::options_description general_desc() {HERE;
    po::options_description description("General");
    description.add_options()
        ("help,h", po::value<bool>()->default_value(false)->implicit_value(true), "print this help")
        ("verbose,v", po::value<bool>()->default_value(false)->implicit_value(true), "verbose output")
        ("test", po::value<int>()->default_value(0)->implicit_value(1));
    return description;
}

po::options_description Program::options_desc() {HERE;
    po::options_description description("Program");
    description.add_options()
      ("popsize,n", po::value(&popsize_)->default_value(popsize_))
      ("initial,q", po::value(&initial_freq_)->default_value(initial_freq_))
      ("generations,g", po::value(&num_generations_)->default_value(num_generations_))
      ("outdir,o", po::value(&outdir_)->default_value(outdir_));
    description.add(Haploid::options_desc());
    description.add(Transposon::options_desc());
    return description;
}

po::options_description Program::positional_desc() {HERE;
    po::options_description description("Positional");
    description.add_options()
        ("nrep", po::value(&num_repeats_)->default_value(num_repeats_));
    return description;
}

void Program::help_and_exit() {HERE;
    auto description = general_desc();
    description.add(options_desc());
    // do not print positional arguments as options
    std::cout << "Usage: tek [options]\n" << std::endl;
    description.print(std::cout);
    throw wtl::ExitSuccess();
}

//! Unit test for each class
inline void test(const int flg) {HERE;
    switch (flg) {
      case 0:
        break;
      case 1:
        Transposon::test();
        Haploid::test();
        Population::test();
        throw wtl::ExitSuccess();
      default:
        throw std::runtime_error("Unknown argument for --test");
    }
}

Program::Program(const std::vector<std::string>& arguments) {HERE;
    std::cout << wtl::join(arguments, " ") << std::endl;
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);

    auto description = general_desc();
    description.add(options_desc());
    description.add(positional_desc());
    po::positional_options_description positional;
    positional.add("nrep", 1);
    po::variables_map vm;
    po::store(po::command_line_parser({arguments.begin() + 1, arguments.end()}).
              options(description).
              positional(positional).run(), vm);
    if (vm["help"].as<bool>()) {help_and_exit();}
    po::notify(vm);
    Transposon::set_parameters();
    Haploid::set_parameters(popsize_, Population::THETA, Population::RHO);

    if (vm["verbose"].as<bool>()) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << wtl::flags_into_string(vm) << std::endl;
    }
    test(vm["test"].as<int>());
}

void Program::run() {HERE;
    if (!outdir_.empty()) wtl::mkdir(outdir_);
    wtl::Pushd cd(outdir_);
    try {
        while (true) {
            Population pop(popsize_, initial_freq_);
            if (pop.evolve(num_generations_)) break;
        }
    } catch (const wtl::KeyboardInterrupt& e) {
        std::cerr << e.what() << std::endl;
    }
}

} // namespace tek
