//
// Created by Jan on 20.04.2022.
//

#include "SimulationEnvironment.h"

#include <iostream>
#include <utility>



SimulationEnvironment::SimulationEnvironment(json initialConfig, InfectionRegime& infectionRegieme) : infectionRegieme(infectionRegieme) {
    config = std::move(initialConfig);
    generation = 0;
}

void SimulationEnvironment::initializeHostAllelePool() {
    hostAllelePool.alleles.resize(config["hosts"]["species_n"]);
    for( int species_i = 0; species_i < config["hosts"]["species_n"]; species_i ++ ){
        hostAllelePool.alleles[species_i].reserve(config["hosts"]["alleles_per_species_initial"]);
        for(int i = 0; i < config["hosts"]["alleles_per_species_initial"]; i++){
            hostAllelePool.alleles[species_i].emplace_back(Allele(i, Helper::gen_random(config["hosts"]["allele_sequence_length"])));
        }
    }
}

void SimulationEnvironment::initializePathogenAllelePool() {
    pathogenAllelePool.alleles.resize(config["pathogens"]["species_n"]);
    for( int species_i = 0; species_i < config["pathogens"]["species_n"]; species_i++ ){
        pathogenAllelePool.alleles[species_i].reserve(config["pathogens"]["haplotypes_per_species_initial"]);
        for(int i = 0; i < config["pathogens"]["haplotypes_per_species_initial"]; i++){
            pathogenAllelePool.alleles[species_i].emplace_back(Allele(i, Helper::gen_random(config["pathogens"]["haplotype_sequence_length"])));
        }
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
    hostPool.hosts.resize(config["hosts"]["species_n"]);
    for( int species_i = 0; species_i < config["hosts"]["species_n"]; species_i ++) {
        hostPool.hosts[species_i].reserve(config["hosts"]["n"]);
        for( int i = 0; i < config["hosts"]["n"]; i++ ){
            hostPool.hosts[species_i].emplace_back(Host(i, 1, species_i));
            for(int j = 0; j < config["hosts"]["genes_per_chromosome_initial"]; j++){
                int randomAlleleId_1 = rand() % hostAllelePool.alleles[species_i].size();
                int randomAlleleId_2 = rand() % hostAllelePool.alleles[species_i].size();
                hostPool.hosts[species_i][i].chromosome_1_allele_ids.emplace_back(randomAlleleId_1);
                hostPool.hosts[species_i][i].chromosome_2_allele_ids.emplace_back(randomAlleleId_2);
            }
        }
    }
}

void SimulationEnvironment::initializePathogenPool() {
    pathogenPool.pathogens.resize(config["pathogens"]["species_n"]);

    for( int species_i = 0; species_i < config["pathogens"]["species_n"]; species_i ++ ){
        pathogenPool.pathogens[species_i].reserve(config["pathogens"]["n"]);
        for( int i = 0; i < config["pathogens"]["n"]; i++ ){
            int randomHaplotypeId = rand() % pathogenAllelePool.alleles[species_i].size();
            pathogenPool.pathogens[species_i].emplace_back(Pathogen(i,1,species_i,randomHaplotypeId));
        }
    }

}

void SimulationEnvironment::printHost(int species, int index){
    Host& host = hostPool.hosts[species][index];
    host.print();

    std::cout << "  chromosome 1 alleles: " << std::endl;
    for(const auto &alleleId : host.chromosome_1_allele_ids){
        Allele& allele = hostAllelePool.alleles[species][alleleId];
        std::cout << "  allele id: " << allele.id << " -> " << allele.sequence << std::endl;
    }

    std::cout << "  chromosome 2 alleles: " << std::endl;
    for(const auto &alleleId : host.chromosome_2_allele_ids){
        Allele& allele = hostAllelePool.alleles[species][alleleId];
        std::cout << "  allele id: " << allele.id << " -> " << allele.sequence << std::endl;
    }
}

void SimulationEnvironment::printPathogen(int species, int index){
    Pathogen& pathogen = pathogenPool.pathogens[species][index];
    pathogen.print();

    Allele& haplotype = pathogenAllelePool.alleles[species][pathogen.haplotypeId];
    std::cout << "  haplotype id: " << haplotype.id << " -> " << haplotype.sequence << std::endl;
}

void SimulationEnvironment::testMethod() {
    infectionRegieme.testMethod();
}

void SimulationEnvironment::initialize() {
    initializeHostAllelePool();
    initializePathogenAllelePool();

    initializeMeritCache();

    initializeHostPool();
    initializePathogenPool();
}

// implement single simulation step (infection, mutation, reproduction)
void SimulationEnvironment::step() {
    infectionRegieme.infect();
    infectionRegieme.testMethod();
    generation += 1;
}
