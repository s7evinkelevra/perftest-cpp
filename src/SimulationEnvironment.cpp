//
// Created by Jan on 20.04.2022.
//

#include "SimulationEnvironment.h"

#include <iostream>
#include <utility>
#include <random>
#include <chrono>


SimulationEnvironment::SimulationEnvironment(json initialConfig, InfectionRegime& infectionRegieme) : infectionRegieme(infectionRegieme) {
    config = std::move(initialConfig);
    generation = 0;

    // init random number generator
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
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
    for( int host_species_i = 0; host_species_i < config["hosts"]["species_n"]; host_species_i++ ){
        for( auto &hostAllele : hostAllelePool.alleles[host_species_i] ){
            for( int patho_species_i = 0; patho_species_i < config["pathogens"]["species_n"]; patho_species_i++ ){
                for( auto &pathogenAllele : pathogenAllelePool.alleles[patho_species_i] ) {
                    int levDistance = Helper::generate_merit(hostAllele.sequence, pathogenAllele.sequence);
                    meritCache.set(host_species_i, hostAllele.id, patho_species_i, pathogenAllele.id, levDistance);
                }
            }
        }
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

// implement single simulation step
// this represents a single host generation with n pathogen generations
void SimulationEnvironment::step() {
    infectionRegieme.infect();
    infectionRegieme.testMethod();
    generation += 1;

    // get distribution of merits of all allele:haplotype combs
    std::unordered_map<int, int> merit_dist = meritCache.getDistribution();
    for(auto& item: merit_dist){
        std::cout << item.first << " : " << item.second << std::endl;
    }

    // i host generations with j pathogen generations each
    // each host generation contains
    //  n pathogen generations
    //      infection
    //      reproduction
    //      mutation
    //  reproduction
    //  mutation
    for(int pathogen_generation = 0; pathogen_generation < config["infection"]["infections_per_generation"]; pathogen_generation++){

        for(int host_species_index = 0; host_species_index < config["hosts"]["species_n"]; host_species_index++){
            unsigned long host_pop_size = hostPool.hosts[host_species_index].size();
            for(int host_index = 0; host_index < host_pop_size; host_index++) {
                Host& selectedHost = hostPool.hosts[host_species_index][host_index];

                for(int patho_species_index = 0; patho_species_index < config["pathogens"]["species_n"]; patho_species_index++){
                    int selectedPathogenIndex = rand() % pathogenPool.pathogens[patho_species_index].size();
                    Pathogen& selectedPathogen =  pathogenPool.pathogens[patho_species_index][selectedPathogenIndex];

                    // determine smallest lev distance of all host alleles to the selected pathogen haplotype
                    int smallest_lev = 99999;
                    for(const int& allele_id : selectedHost.chromosome_1_allele_ids){
                        const int lev_dist = meritCache.get(host_species_index, allele_id, patho_species_index, selectedPathogen.haplotypeId);
                        if(lev_dist < smallest_lev){
                            smallest_lev = lev_dist;
                        }
                    }
                    for(const int& allele_id : selectedHost.chromosome_2_allele_ids){
                        const int lev_dist = meritCache.get(host_species_index, allele_id, patho_species_index, selectedPathogen.haplotypeId);
                        if(lev_dist < smallest_lev){
                            smallest_lev = lev_dist;
                        }
                    }
                    // lev is lte than the provided threshold, at least one peptide of the pathogens haplotype is successfully presented by at least one mhc allele of the given host on either chromosome
                    if(smallest_lev <= config["infection"]["threshold"]){
                        selectedHost.antigen_presentation_count++;
                        selectedPathogen.no_infection_count++;
                    }else{
                        selectedHost.no_antigen_presentation_count++;
                        selectedPathogen.infection_count++;
                    }

                }
            }
        }
        // pathogen
        std::uniform_real_distribution<double> unif(0,1);
        for(int patho_species_index = 0; patho_species_index < pathogenPool.pathogens.size(); patho_species_index++){
            for(Pathogen& currentPathogen : pathogenPool.pathogens[patho_species_index]){

                //std::string newSequence = pathogenAllelePool.alleles[patho_species_index][currentPathogen.haplotypeId].sequence;

                //for(int haplotype_position = 0; haplotype_position < config["pathogens"]["haplotype_sequence_length"]; haplotype_position++){
                    //if(unif(rng) < config["pathogens"]["mutation_rate_per_peptide"]){
                        // mutate
                        //std::cout << "mutation occured" << std::endl;

                        //newSequence[haplotype_position] = config["aminoacids"].get<std::string>();
                        //int newHaplotypeId = addPathogenAllele(patho_species_index);

                    //}
                //}
            }
        }
        // pathogen reproduction

    }

}

int SimulationEnvironment::addHostAllele(int host_species_id, const std::string& sequence) {
    int newAlleleId = hostAllelePool.addAllele(host_species_id, sequence);

    for(int patho_species_id = 0; patho_species_id < config["pathogens"]["species_n"]; patho_species_id++){
        for(Allele& haplotype : pathogenAllelePool.alleles[patho_species_id]){
            int levDist = Helper::generate_merit(sequence, haplotype.sequence);
            meritCache.set(host_species_id, newAlleleId, patho_species_id, haplotype.id, levDist);
        }
    }
    return newAlleleId;
}

int SimulationEnvironment::addPathogenAllele(int patho_species_id, const std::string& sequence) {
    int newAlleleId = pathogenAllelePool.addAllele(patho_species_id, sequence);

    for(int host_species_id = 0; host_species_id < config["pathogens"]["species_n"]; host_species_id++){
        for(Allele& allele : hostAllelePool.alleles[host_species_id]){
            int levDist = Helper::generate_merit(sequence, allele.sequence);
            meritCache.set(host_species_id, newAlleleId, host_species_id, allele.id, levDist);
        }
    }
    return newAlleleId;
}

std::string SimulationEnvironment::generateSequence(int length) {
    std::string AS = config["aminoacids"];
    std::string tmp_s;
    tmp_s.reserve(length);

    for (int i = 0; i < length; ++i) {
        tmp_s += AS[rand() % (AS.length() - 1)];
    }

    return tmp_s;
}