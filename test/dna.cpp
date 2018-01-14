#include "dna.hpp"

#include <random>

int main() {
    std::mt19937 generator;
    constexpr size_t n = 30u;
    tek::DNA x(n);
    tek::DNA y(n);
    std::cerr << x << std::endl;
    for (size_t i=0u; i<n; ++i) {
        y.flip(i, generator);
    }
    std::cerr << y << std::endl;
    std::cerr << wtl::count(x != y) << std::endl;

    tek::Homolog counter(n);
    counter += x;
    counter += y;
    counter += y;
    std::cerr << tek::DNA(counter.majority()) << std::endl;
}
