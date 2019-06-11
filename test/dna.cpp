#include "dna.hpp"

#include <random>

int main() {
    std::mt19937_64 engine;
    constexpr uint_fast32_t n = 30u;
    tek::DNA<n> x;
    tek::DNA<n> y;
    std::cerr << x << std::endl;
    for (uint_fast32_t i=0u; i<n; ++i) {
        y.flip(i, engine);
    }
    std::cerr << y << std::endl;
    std::cerr << (x - y) << std::endl;

    tek::Homolog<n> counter;
    counter += x;
    counter += y;
    counter += y;
    std::cerr << tek::DNA<n>(counter.majority()) << std::endl;
}
