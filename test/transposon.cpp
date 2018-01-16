#include "transposon.hpp"

#include <wtl/iostr.hpp>

#include <random>
#include <iostream>

void activity_function() {
    auto ost = wtl::make_ofs("tek-activity_function.tsv");
    ost << "alpha\tbeta\tidentity\tactivity\n";
    tek::Transposon::write_activity(ost, 0.70,  6);
    tek::Transposon::write_activity(ost, 0.75, 12);
    tek::Transposon::write_activity(ost, 0.80, 24);
    tek::Transposon::write_activity(ost, 0.85, 48);
    /* R
    read_tsv('tek-activity_function.tsv') %>%
    {ggplot(., aes(identity, activity, group=alpha, colour=alpha)) + geom_line()} %>%
    {ggsave('activity_function.pdf', ., width=5, height=3)}
    */
}

int main() {
    std::mt19937 mt(std::random_device{}());
    tek::Transposon::initialize();
    tek::Transposon wt, mut;
    mut.mutate(mt);
    mut.mutate(mt);
    mut.write_fasta(std::cout);
    std::cout << mut << std::endl;
    mut.indel();
    std::cout << mut << std::endl;
    std::cout << "Hamming distance: " << mut - wt << std::endl;
    tek::TransposonFamily family;
    family += wt;
    family += wt;
    family += mut;
    family.majority().write_fasta(std::cout);
    activity_function();
    return 0;
}
