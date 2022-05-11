//
// Created by Jan on 17.04.2022.
//

#ifndef PERFTEST_CPP_HOSTPOOL_H
#define PERFTEST_CPP_HOSTPOOL_H


#include <vector>
#include <unordered_map>
#include "Host.h"

class HostPool {
public:
/*    HostPool(int hostPoolSize);
    virtual ~HostPool();*/
    std::vector<std::vector<Host>> hosts;
    std::vector<std::unordered_map<int, int>> getAlleleDistributions();
    std::unordered_map<int, int> getAlleleDistribution(int speciesId);

    std::vector<int> getAlleleCounts();
    int getAlleleCount(int speciesId);

    std::vector<double> fitness_sum;

    // calculates all the fitness values for each host and saves it to that hosts fitness
    // also updates the fitness sum for that species
    void updateFitness();


};


#endif //PERFTEST_CPP_HOSTPOOL_H
