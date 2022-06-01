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
    // total of all alleles generated for a species. Also includes alleles that have been purged from the allele pool
    std::vector<int> total_allele_counts_per_species;
    std::vector<std::unordered_map<int, Allele>> alleles;
    unsigned long addAllele(int species_id,int parent_id, int created_at_generation, std::string sequence);

    // remove unused alleles from the allele pool of a species
    void purgeUnused(int species_id, std::unordered_map<int, int>& allele_dist);
};


#endif //PERFTEST_CPP_ALLELEPOOL_H
