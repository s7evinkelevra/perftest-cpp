//
// Created by Jan on 18.04.2022.
//

#include "Pathogen.h"
#include <string>
#include <utility>
#include <iostream>

Pathogen::Pathogen(int pathogenId, double initialFitness, int speciesId, std::string pathogenHaplotype) {
    id = pathogenId;
    species = speciesId;
    fitness = initialFitness;
    haplotype = std::move(pathogenHaplotype);
}

void Pathogen::print() {
    std::cout << "Pathogen \nid: " << id << "\nspecies id: " << species << "\nfitness: " << fitness << "\n  haplotype: " << haplotype << std::endl;
}