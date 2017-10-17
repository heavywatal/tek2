// -*- mode: c++; coding: utf-8 -*-
/*! @file program.cpp
    @brief Implementation of Program class
    @defgroup params Parameters
*/
#include "version.hpp"
#include "program.hpp"
#include "population.hpp"
#include "haploid.hpp"
#include "transposon.hpp"
#include "dna.hpp"

#include <sfmt.hpp>
#include <wtl/exception.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/filesystem.hpp>
#include <wtl/getopt.hpp>

#include <iostream>

namespace tek {

namespace po = boost::program_options;

//! options description for general arguments
inline po::options_description general_desc() {HERE;
    po::options_description description("General");
    description.add_options()
        ("help,h", po::bool_switch(), "print this help")
        ("verbose,v", po::bool_switch(), "verbose output")
        ("test", po::value<int>()->default_value(0)->implicit_value(1));
    return description;
}

/*! @ingroup params

    Command line option | Symbol  | Variable
    ------------------- | ------- | -------------------------
    `-n,--popsize`      | \f$N\f$ | Program::popsize_
    `-q,--initial`      |         | Program::initial_freq_
    `-g,--generations`  |         | Program::num_generations_
    `-s,--split`        |         | Program::num_generations_after_split_
    `-i,--interval`     |         | Program::record_interval_
    `-r,--record`       |         | Program::record_flags_
    `-j,--parallel`     |         | Program::concurrency_
    `-o,--outdir`       |         | Program::outdir_
*/
po::options_description Program::options_desc() {HERE;
    po::options_description description("Program");
    description.add_options()
      ("popsize,n", po::value(&popsize_)->default_value(popsize_))
      ("initial,q", po::value(&initial_freq_)->default_value(initial_freq_))
      ("generations,g", po::value(&num_generations_)->default_value(num_generations_))
      ("split,s", po::value(&num_generations_after_split_)->default_value(num_generations_after_split_))
      ("interval,i", po::value(&record_interval_)->default_value(record_interval_))
      ("record,r", po::value(&record_flags_)->default_value(record_flags_))
      ("parallel,j", po::value(&concurrency_)->default_value(concurrency_))
      ("outdir,o", po::value(&outdir_)->default_value(outdir_));
    description.add(Haploid::options_desc());
    description.add(Transposon::options_desc());
    return description;
}

void Program::help_and_exit() {HERE;
    auto description = general_desc();
    description.add(options_desc());
    // do not print positional arguments as options
    std::cout << "commit " << GIT_COMMIT_HASH
              << " [" << GIT_BRANCH << "]\n"
              << "Date:  " << GIT_COMMIT_TIME << std::endl;
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
        DNA_test(wtl::sfmt());
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
    po::variables_map vm;
    po::store(po::command_line_parser({arguments.begin() + 1, arguments.end()}).
              options(description).run(), vm);
    if (vm["help"].as<bool>()) {help_and_exit();}
    po::notify(vm);
    Transposon::set_parameters();
    Haploid::set_parameters(popsize_, Population::THETA, Population::RHO);

    config_string_ = wtl::flags_into_string(vm);
    if (vm["verbose"].as<bool>()) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << config_string_ << std::endl;
    }
    test(vm["test"].as<int>());
}

void Program::run() {HERE;
    try {
        main();
    } catch (const wtl::KeyboardInterrupt& e) {
        std::cerr << e.what() << std::endl;
    }
}

void Program::main() {HERE;
    wtl::ChDir cd_outdir(outdir_, true);
    while (true) {
        Population pop(popsize_, initial_freq_, concurrency_);
        auto flags = static_cast<Recording>(record_flags_);
        bool good = pop.evolve(num_generations_, record_interval_, flags);
        if (!good) continue;
        wtl::opfstream{"program_options.conf"} << config_string_;
        {
          wtl::ozfstream ost("summary.json.gz");
          pop.write_summary(ost);
        }
        {
          wtl::ozfstream ost("sequence.fa.gz");
          pop.write_fasta(ost);
        }
        if (num_generations_after_split_ == 0U) break;
        Population pop2(pop);
        {
          wtl::ChDir cd("population_1", true);
          good = pop.evolve(num_generations_after_split_, record_interval_, Recording::sequence);
        }
        if (!good) continue;
        {
          wtl::ChDir cd("population_2", true);
          good = pop2.evolve(num_generations_after_split_, record_interval_, Recording::sequence);
        }
        if (!good) continue;
        break;
    }
}

} // namespace tek
