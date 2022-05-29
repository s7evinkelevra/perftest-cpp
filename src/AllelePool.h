//
// Created by Jan on 18.04.2022.
//

#ifndef PERFTEST_CPP_ALLELEPOOL_H
#define PERFTEST_CPP_ALLELEPOOL_H

#include <vector>
#include <unordered_map>
#include "Allele.h"

class AllelePool {
public:
    std::vector<int> counts;
    std::vector<std::vector<Allele>> alleles;
    unsigned long addAllele(int species_id,int parent_id, std::string sequence);

    void purgeUnused(int species_id, std::unordered_map<int, int>& allele_dist);
};


#endif //PERFTEST_CPP_ALLELEPOOL_H
