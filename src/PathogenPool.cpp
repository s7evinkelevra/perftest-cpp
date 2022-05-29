//
// Created by Jan on 18.04.2022.
//

#include "PathogenPool.h"

void PathogenPool::updateFitness() {
    for(int species_i = 0; species_i < pathogens.size(); species_i++){
        double current_fitness_sum = 0;
        for(auto& pathogen : pathogens[species_i]){
            pathogen.updateFitness();
            current_fitness_sum += pathogen.fitness;
        }
        fitness_sum[species_i] = current_fitness_sum;
    }
}

std::unordered_map<int, int> PathogenPool::getHaplotypeDistribution(int speciesId) {
    std::unordered_map<int, int> dist;
    for(auto& pathogen : pathogens[speciesId]){
        dist[pathogen.haplotype_id]++;
    }
    return dist;
}

std::vector<std::unordered_map<int, int>> PathogenPool::getHaplotypeDistributionsPerSpecies() {
    std::vector<std::unordered_map<int,int>> dist_vec;

    for(int species_id = 0; species_id < pathogens.size(); species_id++){
        std::unordered_map<int, int> species_dist = getHaplotypeDistribution(species_id);
        dist_vec.push_back(species_dist);
    }

    return dist_vec;
}
