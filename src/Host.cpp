//
// Created by Jan on 18.04.2022.
//

#include <iostream>
#include "Host.h"

Host::Host(int hostId, double initialFitness, int speciesId) {
    id = hostId;
    fitness = initialFitness;
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
