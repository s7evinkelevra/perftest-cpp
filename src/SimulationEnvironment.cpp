//
// Created by Jan on 20.04.2022.
//

#include "SimulationEnvironment.h"

#include <iostream>
#include <utility>



SimulationEnvironment::SimulationEnvironment(json initialConfig, InfectionRegieme& infectionRegieme) : infectionRegieme(infectionRegieme) {
    config = std::move(initialConfig);
}

void SimulationEnvironment::initializeHostAllelePool() {
    hostAllelePool.alleles.reserve(config["hosts"]["alleles_per_species_initial"]);
    for(int i = 0; i < config["hosts"]["alleles_per_species_initial"]; i++){
        hostAllelePool.alleles.emplace_back(Allele(i, Helper::gen_random(config["hosts"]["allele_sequence_length"])));
    }
}

void SimulationEnvironment::initializePathogenAllelePool() {
    pathogenAllelePool.alleles.reserve(config["pathogens"]["haplotypes_per_species_initial"]);
    for(int i = 0; i < config["pathogens"]["haplotypes_per_species_initial"]; i++){
        pathogenAllelePool.alleles.emplace_back(Allele(i, Helper::gen_random(config["pathogens"]["haplotype_sequence_length"])));
    }
}

void SimulationEnvironment::initializeMeritCache() {
    for( auto &hostAllele : hostAllelePool.alleles ){
        std::deque<int> row;
        for( auto &pathogenAllele : pathogenAllelePool.alleles ) {
            int levDistance = Helper::generate_merit(hostAllele.sequence, pathogenAllele.sequence);
            row.push_back(levDistance);
        }
        meritCache.cache.push_back(row);
    }
};

void SimulationEnvironment::initializeHostPool() {
    hostPool.hosts.reserve(config["hosts"]["n"]);
    for(int i = 0; i < config["hosts"]["n"]; i++){
        hostPool.hosts.emplace_back(Host(i, 1, 0));
        for(int j = 0; j < config["hosts"]["genes_per_chromosome_initial"]; j++){
            int randomAlleleId_1 = rand() % hostAllelePool.alleles.size();
            int randomAlleleId_2 = rand() % hostAllelePool.alleles.size();
            hostPool.hosts[i].chromosome_1_allele_ids.emplace_back(randomAlleleId_1);
            hostPool.hosts[i].chromosome_2_allele_ids.emplace_back(randomAlleleId_2);
        }
    }
}

void SimulationEnvironment::initializePathogenPool() {
    pathogenPool.pathogens.reserve(config["pathogens"]["n"]);
    for(int i = 0; i < config["pathogens"]["n"]; i++){
        int randomHaplotypeId = rand() % pathogenAllelePool.alleles.size();
        pathogenPool.pathogens.emplace_back(Pathogen(i,1,0,randomHaplotypeId));
    }

}

void SimulationEnvironment::printHost(int index){
    Host& host = hostPool.hosts[index];
    host.print();

    std::cout << "  chromosome 1 alleles: " << std::endl;
    for(const auto &alleleId : host.chromosome_1_allele_ids){
        Allele& allele = hostAllelePool.alleles[alleleId];
        std::cout << "  allele id: " << allele.id << " -> " << allele.sequence << std::endl;
    }

    std::cout << "  chromosome 2 alleles: " << std::endl;
    for(const auto &alleleId : host.chromosome_2_allele_ids){
        Allele& allele = hostAllelePool.alleles[alleleId];
        std::cout << "  allele id: " << allele.id << " -> " << allele.sequence << std::endl;
    }
}

void SimulationEnvironment::printPathogen(int index){
    Pathogen& pathogen = pathogenPool.pathogens[index];
    pathogen.print();

    Allele& haplotype = pathogenAllelePool.alleles[pathogen.haplotypeId];
    std::cout << "  haplotype id: " << haplotype.id << " -> " << haplotype.sequence << std::endl;
}

