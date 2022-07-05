//
// Created by Jan on 18.04.2022.
//

#include <iostream>
#include "Host.h"

Host::Host(int parentId1, double parentFitness1, int parentId2, double parentFitness2, int hostId, double fitnessMinimum, int speciesId) {
    parent_id_1 = parentId1;
    parent_fitness_1 = parentFitness1;
    parent_id_2 = parentId2;
    parent_fitness_2 = parentFitness2;
    id = hostId;
    fitness = fitnessMinimum;
    fitness_minimum = fitnessMinimum;
    species = speciesId;
    chromosome_1_allele_ids.reserve(20);
    chromosome_2_allele_ids.reserve(20);

    antigen_presentation_count = 0;
    no_antigen_presentation_count = 0;
    //std::cout << "initializing host "<< id << " \n";
}

void Host::print() {
    std::cout << "Host\nid: " << id << "\nspecies id: " << species << "\nfitness: " << fitness << "\nsuccessfull antigen presentations: " << antigen_presentation_count << "\nno antigen presentations: " << no_antigen_presentation_count << std::endl;
/*  std::cout << "\n Chromosome 1 allele Ids: ";
    for(const auto &alleleId : chromosome_1_allele_ids){
        std::cout << alleleId << " ";
    }
    std::cout << "\n Chromosome 2 allele Ids: ";
    for(const auto &alleleId : chromosome_2_allele_ids){
        std::cout << alleleId << " ";
    }
    std::cout << std::endl;*/
}

void Host::updateFitness() {
    if(antigen_presentation_count + no_antigen_presentation_count == 0){
        fitness = fitness_minimum;
    }else{
        double new_fitness = (double)antigen_presentation_count / ((double)antigen_presentation_count + (double)no_antigen_presentation_count);
        fitness = std::max(fitness_minimum, new_fitness);
    }
}

void Host::resetAntigenPresentationCount(int initialFitness) {
    fitness = initialFitness;
    antigen_presentation_count = 0;
    no_antigen_presentation_count = 0;
}
