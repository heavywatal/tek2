/*! @file program.cpp
    @brief Implementation of Program class
    @defgroup params Parameters
*/
#include "version.hpp"
#include "program.hpp"
#include "population.hpp"
#include "haploid.hpp"
#include "transposon.hpp"

#include <wtl/exception.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/chrono.hpp>
#include <wtl/zlib.hpp>
#include <wtl/filesystem.hpp>
#include <clippson/clippson.hpp>

namespace tek {

//! variables map
nlohmann::json VM;

//! Options description for general purpose
inline clipp::group general_options(nlohmann::json* vm) {HERE;
    return (
      wtl::option(vm, {"h", "help"}, false, "print this help"),
      wtl::option(vm, {"version"}, false, "print version"),
      wtl::option(vm, {"v", "verbose"}, false, "verbose output")
    ).doc("General:");
}

/*! @ingroup params

    Command line option | Symbol  | Variable
    ------------------- | ------- | -------------------------
    `-n,--popsize`      | \f$N\f$ |
    `-q,--initial`      |         |
    `-g,--generations`  |         |
    `-s,--split`        |         |
    `-i,--interval`     |         |
    `-r,--record`       |         |
    `-o,--outdir`       |         |
*/
inline clipp::group program_options(nlohmann::json* vm) {HERE;
    const std::string outdir = wtl::strftime("tek_%Y%m%d_%H%M%S");
    const int seed = static_cast<int>(std::random_device{}()); // 32-bit signed integer for R
    return (
      wtl::option(vm, {"n", "popsize"}, 500u),
      wtl::option(vm, {"q", "initial"}, 1u,
        "initial number of individuals with TE"),
      wtl::option(vm, {"g", "generations"}, 1000u,
        "maximum number of generations to simulate"),
      wtl::option(vm, {"H", "hyperactivate"}, 0u,
        "time to introduce a hyperactivating mutation"),
      wtl::option(vm, {"s", "split"}, 0u,
        "number of generations to simulate after population split"),
      wtl::option(vm, {"i", "interval"}, 10u,
        "interval of recording"),
      wtl::option(vm, {"r", "record"}, 3,
        "enum Recording"),
      wtl::option(vm, {"o", "outdir"}, outdir),
      wtl::option(vm, {"seed"}, seed)
    ).doc("Program:");
}

/*! @ingroup params

    Command line option | Symbol        | Variable
    ------------------- | ------------- | -------------------------
    `--sample`          |               | PopulationParams::SAMPLE_SIZE
    `-j,--parallel`     |               | PopulationParams::CONCURRENCY
    `-c,--coexist`      |               | PopulationParams::MAX_COEXISTENCE
*/
inline clipp::group
population_options(nlohmann::json* vm, PopulationParams* p) {HERE;
    return (
      wtl::option(vm, {"sample"}, &p->SAMPLE_SIZE),
      wtl::option(vm, {"j", "parallel"}, &p->CONCURRENCY),
      wtl::option(vm, {"c", "coexist"}, &p->MAX_COEXISTENCE)
    ).doc("Population:");
}

/*! @ingroup params

    Command line option | Symbol        | Variable
    ------------------- | ------------- | -------------------------
    `--xi`              | \f$\xi\f$     | HaploidParams::XI
    `--nu`              | \f$\nu\f$     | HaploidParams::EXCISION_RATE
    `--lambda`          | \f$\lambda\f$ | HaploidParams::MEAN_SELECTION_COEF
*/
inline clipp::group
haploid_options(nlohmann::json* vm, HaploidParams* p) {HERE;
    return (
      wtl::option(vm, {"xi"}, &p->XI),
      wtl::option(vm, {"nu"}, &p->EXCISION_RATE),
      wtl::option(vm, {"lambda"}, &p->MEAN_SELECTION_COEF)
    ).doc("Haploid:");
}

/*! @ingroup params

    Command line option | Symbol        | Variable
    ------------------- | ------------- | -------------------------
    `-a,--alpha`        | \f$\alpha\f$  | TransposonParams::ALPHA
    `-b,--beta`         | \f$\beta\f$   | TransposonParams::BETA
    `-l,--lower`        | \f$d_l\f$     | TransposonParams::LOWER_THRESHOLD
    `-u,--upper`        | \f$d_u\f$     | TransposonParams::UPPER_THRESHOLD
*/
inline clipp::group
transposon_options(nlohmann::json* vm, TransposonParams* p) {HERE;
    return (
      wtl::option(vm, {"a", "alpha"}, &p->ALPHA),
      wtl::option(vm, {"b", "beta"}, &p->BETA),
      wtl::option(vm, {"l", "lower"}, &p->LOWER_THRESHOLD),
      wtl::option(vm, {"u", "upper"}, &p->UPPER_THRESHOLD)
    ).doc("Transposon:");
}

Program::Program(const std::vector<std::string>& arguments) {HERE;
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);

    PopulationParams population_params;
    HaploidParams haploid_params;
    TransposonParams transposon_params;

    VM.clear();
    nlohmann::json vm_local;
    auto cli = (
      general_options(&vm_local),
      program_options(&VM),
      population_options(&VM, &population_params),
      haploid_options(&VM, &haploid_params),
      transposon_options(&VM, &transposon_params)
    );
    wtl::parse(cli, arguments);
    auto fmt = wtl::doc_format();
    if (vm_local["help"]) {
        std::cout << "Usage: " << PROJECT_NAME << " [options]\n\n";
        std::cout << clipp::documentation(cli, fmt) << "\n";
        throw wtl::ExitSuccess();
    }
    if (vm_local["version"]) {
        std::cout << PROJECT_VERSION << "\n";
        throw wtl::ExitSuccess();
    }
    Population::param(population_params);
    Haploid::param(haploid_params);
    Transposon::param(transposon_params);
    config_ = VM.dump(2);
    if (vm_local["verbose"]) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << config_ << std::endl;
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
    const size_t popsize_ = VM.at("popsize");
    const size_t initial_freq_ = VM.at("initial");
    const size_t num_generations_ = VM.at("generations");
    const size_t hyperactivate = VM.at("hyperactivate");
    const size_t num_generations_after_split_ = VM.at("split");
    const size_t record_interval_ = VM.at("interval");
    const int record_flags_ = VM.at("record");
    const std::string outdir_ = VM.at("outdir");
    Population::seed(VM.at("seed"));
    wtl::ChDir cd_outdir(outdir_, true);
    while (true) {
        Population pop(popsize_, initial_freq_);
        auto flags = static_cast<Recording>(record_flags_);
        bool good = pop.evolve(num_generations_, record_interval_, flags, hyperactivate);
        if (!good) continue;
        wtl::make_ofs("config.json") << config_;
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
