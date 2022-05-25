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

std::vector<std::unordered_map<int, int>> HostPool::getAlleleDistributions(int locus_id) {
    std::vector<std::unordered_map<int,int>> dist_vec;

    for(int species_id = 0; species_id < hosts.size(); species_id++){
        std::unordered_map<int, int> species_dist = getAlleleDistribution(species_id, locus_id);
        dist_vec.push_back(species_dist);
    }

    return dist_vec;
}

std::unordered_map<int, int> HostPool::getAlleleDistribution(int speciesId, int locus_id) {
    std::unordered_map<int, int> dist;

    for(auto& host : hosts[speciesId]){
        if(host.chromosome_1_allele_ids.size() > locus_id){
            dist[host.chromosome_1_allele_ids[locus_id]]++;
        };

        if(host.chromosome_2_allele_ids.size() > locus_id){
            dist[host.chromosome_2_allele_ids[locus_id]]++;
        };
    }

    return dist;
}

std::vector<int> HostPool::getAlleleCounts(int locus_id) {
    std::vector<int> counts;
    counts.reserve(hosts.size());
    for(int species_id = 0; species_id < hosts.size(); species_id++){
        counts.push_back(getAlleleCount(species_id, locus_id));
    }
    return counts;
}

int HostPool::getAlleleCount(int speciesId, int locus_id) {
    return (int) getAlleleDistribution(speciesId, locus_id).size();
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
