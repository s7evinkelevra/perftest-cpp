//
// Created by Jan on 18.04.2022.
//

#ifndef PERFTEST_CPP_PATHOGENPOOL_H
#define PERFTEST_CPP_PATHOGENPOOL_H

#include <vector>
#include <unordered_map>
#include "Pathogen.h"

class PathogenPool {
public:
    std::vector<std::vector<Pathogen>> pathogens;

    // initialized in the simulation environment
    std::vector<double> fitness_sum;
    void updateFitness();

    std::unordered_map<int, int> getHaplotypeDistribution(int speciesId);
    std::vector<std::unordered_map<int, int>> getHaplotypeDistributions();


};


#endif //PERFTEST_CPP_PATHOGENPOOL_H
