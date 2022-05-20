//
// Created by Jan on 18.04.2022.
//

#include "Pathogen.h"
#include <string>
#include <utility>
#include <iostream>

Pathogen::Pathogen(int parentId, int pathogenId, double initialFitness, int speciesId, int pathogenHaplotypeId) {
    parent_id = parentId;
    id = pathogenId;
    species = speciesId;
    fitness = initialFitness;
    haplotype_id = pathogenHaplotypeId;

    infection_count = 0;
    no_infection_count = 0;
}

void Pathogen::print() {
    std::cout << "Pathogen \nid: " << id << "\nparent id:" << parent_id << "\nspecies id: " << species << "\nfitness: " << fitness << "\ninfection count: " << infection_count << "\nno infection count: " << no_infection_count << std::endl;
}

void Pathogen::updateFitness() {
    if(infection_count + no_infection_count == 0){
        fitness = 0;
    }else{
        fitness = (double)infection_count / ((double)infection_count + (double)no_infection_count);
    }
}

void Pathogen::resetInfectionCount(int initialFitness) {
    fitness = initialFitness;
    infection_count = 0;
    no_infection_count = 0;
}
