//
// Created by Jan on 18.04.2022.
//

#include "Pathogen.h"
#include <string>
#include <utility>
#include <iostream>

Pathogen::Pathogen(int pathogenId, double initialFitness, int speciesId, int pathogenHaplotypeId) {
    id = pathogenId;
    species = speciesId;
    fitness = initialFitness;
    haplotypeId = pathogenHaplotypeId;

    infection_count = 0;
    no_infection_count = 0;
}

void Pathogen::print() {
    std::cout << "Pathogen \nid: " << id << "\nspecies id: " << species << "\nfitness: " << fitness << "\ninfection count: " << infection_count << "\nno infection count: " << no_infection_count << std::endl;
}

void Pathogen::updateFitness() {
    if(infection_count + no_infection_count == 0){
        fitness = 0;
    }else{
        fitness = (double)infection_count / ((double)infection_count + (double)no_infection_count);
    }
}
