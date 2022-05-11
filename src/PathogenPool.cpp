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
