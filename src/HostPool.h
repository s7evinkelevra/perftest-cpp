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

    // return allele counts per locus
    // e.g. host[0].chromosome_1 -> alleles [2,5]
    // e.g. host[0].chromosome_2 -> alleles [5,5]
    // dist for locus 0
    // dist[2] = 1
    // dist[5] = 1
    // dist for locus 1
    // dist[2] = 0
    // dist[5] = 2
    std::vector<std::vector<Host>> hosts;
    std::vector<std::vector<std::unordered_map<int, int>>> getAlleleDistributionsPerSpeciesPerLocus();
    std::vector<std::unordered_map<int, int>> getAlleleDistributionsPerLocus(int species_id);
    std::unordered_map<int, int> getAlleleDistribution(int species_id, int locus_id);

    // returns allele counts ACROSS ALL LOCI
    // e.g. host[0].chromosome_1 -> alleles [2,5]
    // e.g. host[0].chromosome_2 -> alleles [5,5]
    // distribution for both loci:
    // dist[2] = 1
    // dist[5] = 3
    std::vector<std::unordered_map<int, int>> getAlleleDistributionAcrossAllLociPerSpecies();
    std::unordered_map<int, int> getAlleleDistributionAcrossAllLoci(int species_id);

    std::vector<int> getAlleleCounts(int locus_id);
    int getAlleleCount(int speciesId, int locus_id);

    // vector of the sums of fitness for each species
    // updated in updateFitness
    std::vector<double> fitness_sum;

    // calculates all the fitness values for each host and saves it to that hosts fitness
    // also updates the fitness sum for that species
    void updateFitness();

    // manually track the maximum number of loci found in an individual
    int max_loci_count;
    // check all individuals and set max loci count
    void updateMaxLociCount();

};


#endif //PERFTEST_CPP_HOSTPOOL_H
