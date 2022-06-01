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

std::vector<std::vector<std::unordered_map<int, int>>> HostPool::getAlleleDistributionsPerSpeciesPerLocus() {
    std::vector<std::vector<std::unordered_map<int,int>>> species_loci_dist_vec;
    species_loci_dist_vec.reserve(hosts.size());

    for(int species_i = 0; species_i < hosts.size(); species_i++){
        std::vector<std::unordered_map<int,int>> species_dist_vec = getAlleleDistributionsPerLocus(species_i);
        species_loci_dist_vec.emplace_back(species_dist_vec);
    }

    return species_loci_dist_vec;
}


std::vector<std::unordered_map<int, int>> HostPool::getAlleleDistributionsPerLocus(int species_id) {
    std::vector<std::unordered_map<int,int>> locus_dist_vec;
    locus_dist_vec.reserve(max_loci_count);

    for(int locus_i = 0; locus_i < max_loci_count; locus_i++){
        std::unordered_map<int, int> locus_dist = getAlleleDistribution(species_id, locus_i);
        locus_dist_vec.emplace_back(locus_dist);
    }

    return locus_dist_vec;
}

std::unordered_map<int, int> HostPool::getAlleleDistribution(int species_id, int locus_id) {
    std::unordered_map<int, int> dist;

    for(auto& host : hosts[species_id]){
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
        counts.emplace_back(getAlleleCount(species_id, locus_id));
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

void HostPool::updateMaxLociCount() {
    for(int species_i = 0; species_i < hosts.size(); species_i++){
        for(auto& host : hosts[species_i]){
            //max_loci_count = std::max({(int)max_loci_count, (int)host.chromosome_1_allele_ids.size(), (int)host.chromosome_2_allele_ids.size()});
        }
    }
}


std::vector<std::unordered_map<int, int>> HostPool::getAlleleDistributionAcrossAllLociPerSpecies() {
    std::vector<std::unordered_map<int,int>> species_dist_vec;
    species_dist_vec.reserve(hosts.size());

    for(int species_i = 0; species_i < hosts.size(); species_i++){
        std::unordered_map<int, int> dist = getAlleleDistributionAcrossAllLoci(species_i);
        species_dist_vec.emplace_back(dist);
    }

    return species_dist_vec;
}


std::unordered_map<int, int> HostPool::getAlleleDistributionAcrossAllLoci(int species_id) {
    std::unordered_map<int, int> dist;

    for(auto& host : hosts[species_id]){
        for(int& id : host.chromosome_1_allele_ids){
            dist[id]++;
        }
        for(int& id : host.chromosome_2_allele_ids){
            dist[id]++;
        }
    }

    return dist;
}


