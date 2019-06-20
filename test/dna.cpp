#include "dna.hpp"

#include <random>

int main() {
    tek::DNA<4> letters(std::valarray<uint_fast8_t>{0, 1, 2, 3});
    std::cerr << letters << std::endl;

    std::mt19937_64 engine;
    constexpr uint_fast32_t n = 30u;
    tek::DNA<n> x;
    tek::DNA<n> y;
    std::cerr << x << std::endl;
    for (uint_fast32_t i=0u; i<n; ++i) {
        y.flip(i, engine);
    }
    tek::DNA<n> z = y;
    for (uint_fast32_t i=0u; i<n; ++i) {
        z.flip(i, engine);
    }
    std::cerr << y << std::endl;
    std::cerr << z << std::endl;
    std::cerr << (x - y) << std::endl;
    std::cerr << z.count() << std::endl;

    tek::Homolog<n> counter;
    counter += x;
    counter += y;
    counter += z;
    counter += z;
    std::cerr << tek::DNA<n>(counter.majority()) << std::endl;
}
