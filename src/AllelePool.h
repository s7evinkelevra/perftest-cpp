//
// Created by Jan on 18.04.2022.
//

#ifndef PERFTEST_CPP_ALLELEPOOL_H
#define PERFTEST_CPP_ALLELEPOOL_H

#include <vector>
#include "Allele.h"

class AllelePool {
public:
    std::vector<std::vector<Allele>> alleles;
    unsigned long addAllele(int species_id, std::string sequence);
};


#endif //PERFTEST_CPP_ALLELEPOOL_H
