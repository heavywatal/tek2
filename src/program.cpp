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

#include <wtl/exception.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/chrono.hpp>
#include <wtl/zlib.hpp>
#include <wtl/filesystem.hpp>
#include <wtl/getopt.hpp>
#include <boost/program_options.hpp>

namespace tek {

namespace po = boost::program_options;

//! directory for output
const std::string OUT_DIR = wtl::strftime("tek_%Y%m%d_%H%M%S");

//! options description for general arguments
inline po::options_description general_desc() {HERE;
    po::options_description description("General");
    description.add_options()
        ("help,h", po::bool_switch(), "print this help")
        ("verbose,v", po::bool_switch(), "verbose output")
    ;
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
    `-o,--outdir`       |         | Program::outdir_
*/
po::options_description Program::options_desc() {HERE;
    auto po_value = [](auto* x) {return po::value(x)->default_value(*x);};
    po::options_description description("Program");
    description.add_options()
      ("popsize,n", po_value(&popsize_))
      ("initial,q", po_value(&initial_freq_))
      ("generations,g", po_value(&num_generations_))
      ("split,s", po_value(&num_generations_after_split_))
      ("interval,i", po_value(&record_interval_))
      ("record,r", po_value(&record_flags_))
      ("outdir,o", po::value(&outdir_)->default_value(OUT_DIR))
    ;
    description.add(Population::options_desc());
    description.add(Haploid::options_desc());
    description.add(Transposon::options_desc());
    return description;
}

[[noreturn]] void Program::help_and_exit() {HERE;
    auto description = general_desc();
    description.add(options_desc());
    // do not print positional arguments as options
    std::cout << "commit " << GIT_COMMIT_HASH
              << " [" << GIT_BRANCH << "]\n"
              << "Date:  " << GIT_COMMIT_TIME << "\n";
    std::cout << "Usage: tek [options]\n\n";
    description.print(std::cout);
    throw wtl::ExitSuccess();
}

Program::Program(const std::vector<std::string>& arguments) {HERE;
    wtl::join(arguments, std::cout, " ") << std::endl;
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

    config_string_ = wtl::flags_into_string(vm);
    if (vm["verbose"].as<bool>()) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << config_string_ << std::endl;
    }
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
        Population pop(popsize_, initial_freq_);
        auto flags = static_cast<Recording>(record_flags_);
        bool good = pop.evolve(num_generations_, record_interval_, flags);
        if (!good) continue;
        wtl::make_ofs("program_options.conf") << config_string_;
        if (static_cast<bool>(flags & Recording::sequence)) {
            wtl::zlib::ofstream ost("sequence.fa.gz");
            pop.write_fasta(ost);
        }
        if (static_cast<bool>(flags & Recording::summary)) {
            wtl::zlib::ofstream ost("summary.json.gz");
            pop.write_summary(ost);
        }
        if (num_generations_after_split_ == 0u) break;
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
