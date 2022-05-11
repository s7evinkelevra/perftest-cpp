//
// Created by Jan on 18.04.2022.
//

#ifndef PERFTEST_CPP_HOST_H
#define PERFTEST_CPP_HOST_H

#include <vector>
#include "Allele.h"

class Host {
public:
    Host(int hostId, double initialFitness, int speciesId);
    int id;
    double fitness;
    int species;
    std::vector<int> chromosome_1_allele_ids;
    std::vector<int> chromosome_2_allele_ids;

    int antigen_presentation_count;
    int no_antigen_presentation_count;

    void updateFitness();

    void print();
};


#endif //PERFTEST_CPP_HOST_H
