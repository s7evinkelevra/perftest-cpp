//
// Created by Jan on 17.04.2022.
//
#include <iostream>
#include "HostPool.h"
#include "Host.h"

/*
// host pool gets initialized by the simulation environment, not here
HostPool::HostPool() {
    fitness_sum.resize()
}
*/

std::vector<std::unordered_map<int, int>> HostPool::getAlleleDistributions() {
    std::vector<std::unordered_map<int,int>> dist_vec;

    for(int species_id = 0; species_id < hosts.size(); species_id++){
        std::unordered_map<int, int> species_dist = getAlleleDistribution(species_id);
        dist_vec.push_back(species_dist);
    }

    return dist_vec;
}

std::unordered_map<int, int> HostPool::getAlleleDistribution(int speciesId) {
    std::unordered_map<int, int> dist;

    for(auto& host : hosts[speciesId]){
        for(auto& allele_id : host.chromosome_1_allele_ids){
            dist[allele_id]++;
        }
        for(auto& allele_id : host.chromosome_2_allele_ids){
            dist[allele_id]++;
        }
    }

    return dist;
}

std::vector<int> HostPool::getAlleleCounts() {
    std::vector<int> counts;
    counts.reserve(hosts.size());
    for(int species_id = 0; species_id < hosts.size(); species_id++){
        counts.push_back(getAlleleCount(species_id));
    }
    return counts;
}

int HostPool::getAlleleCount(int speciesId) {
    return (int) getAlleleDistribution(speciesId).size();
}

void HostPool::updateFitness() {
    for(int species_i = 0; species_i < hosts.size(); species_i++){
        double current_fitness_sum = 0;
        for(auto& host : hosts[species_i]){
            host.updateFitness();
            current_fitness_sum += host.fitness;
        }
        fitness_sum[species_i] = current_fitness_sum;
    }
}
