#include "population.hpp"

#include <iostream>

int main() {
    tek::Population pop(6, 6);
    std::cout << pop << std::endl;
    pop.evolve(1u, -1u);
    std::cout << pop << std::endl;
    pop.write_summary(std::cout);
    pop.write_fasta(std::cout);
}
