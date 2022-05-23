//
// Created by Jan on 20.04.2022.
//

#include "SimulationEnvironment.h"

#include <iostream>
#include <utility>



SimulationEnvironment::SimulationEnvironment(json initialConfig, InfectionRegime& infectionRegieme) : infectionRegieme(infectionRegieme) {
    config = std::move(initialConfig);
    totalHostGenerations = 0;
    totalPathogenGenerations = 0;

    bInfection = true;
    bHostFitnessproportionalReproduction = true;
    bPathogenFitnessproportionalReproduction = true;
    bHostMutation = true;
    bPathogenMutation = true;

    rng = Random();
}

void SimulationEnvironment::initializeHostAllelePool() {
    hostAllelePool.alleles.resize(config["hosts"]["species_n"]);
    for( int species_i = 0; species_i < config["hosts"]["species_n"]; species_i ++ ){
        hostAllelePool.alleles[species_i].reserve(config["hosts"]["alleles_per_species_initial"]);
        for(int i = 0; i < config["hosts"]["alleles_per_species_initial"]; i++){
            hostAllelePool.alleles[species_i].emplace_back(Allele(i, generateSequence(config["hosts"]["allele_sequence_length"])));
        }
    }
}

void SimulationEnvironment::initializePathogenAllelePool() {
    pathogenAllelePool.alleles.resize(config["pathogens"]["species_n"]);
    for( int species_i = 0; species_i < config["pathogens"]["species_n"]; species_i++ ){
        pathogenAllelePool.alleles[species_i].reserve(config["pathogens"]["haplotypes_per_species_initial"]);
        for(int i = 0; i < config["pathogens"]["haplotypes_per_species_initial"]; i++){
            pathogenAllelePool.alleles[species_i].emplace_back(Allele(i, generateSequence(config["pathogens"]["haplotype_sequence_length"])));
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
    hostPool.fitness_sum.resize(config["hosts"]["species_n"]);
    int initialFitness = config["hosts"]["initial_fitness"];

    for( int species_i = 0; species_i < config["hosts"]["species_n"]; species_i ++) {
        hostPool.hosts[species_i].reserve(config["hosts"]["n"]);
        for( int i = 0; i < config["hosts"]["n"]; i++ ){
            hostPool.hosts[species_i].emplace_back(Host(0, 0, i, initialFitness, species_i));
            for(int j = 0; j < config["hosts"]["genes_per_chromosome_initial"]; j++){
                int randomAlleleId_1 = rng.sampleIntUniUnsignedInt(0,hostAllelePool.alleles[species_i].size() - 1);
                int randomAlleleId_2 = rng.sampleIntUniUnsignedInt(0,hostAllelePool.alleles[species_i].size() - 1);
                hostPool.hosts[species_i][i].chromosome_1_allele_ids.emplace_back(randomAlleleId_1);
                hostPool.hosts[species_i][i].chromosome_2_allele_ids.emplace_back(randomAlleleId_2);
            }
        }
    }
}

void SimulationEnvironment::initializePathogenPool() {
    pathogenPool.pathogens.resize(config["pathogens"]["species_n"]);
    pathogenPool.fitness_sum.resize(config["pathogens"]["species_n"]);
    int initialFitness = config["pathogens"]["initial_fitness"];

    for( int species_i = 0; species_i < config["pathogens"]["species_n"]; species_i ++ ){
        pathogenPool.pathogens[species_i].reserve(config["pathogens"]["n"]);
        for( int i = 0; i < config["pathogens"]["n"]; i++ ){
            int randomHaplotypeId = rng.sampleIntUniUnsignedInt(0, pathogenAllelePool.alleles[species_i].size() - 1);
            pathogenPool.pathogens[species_i].emplace_back(Pathogen(0,i,initialFitness,species_i,randomHaplotypeId));
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

    Allele& haplotype = pathogenAllelePool.alleles[species][pathogen.haplotype_id];
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
// this represents a single host totalHostGenerations with n pathogen generations
void SimulationEnvironment::step() {
    infectionRegieme.infect();
    infectionRegieme.testMethod();
    totalHostGenerations++;

    // get distribution of merits of all allele:haplotype combs
//    std::unordered_map<int, int> merit_dist = meritCache.getDistribution();
//    std::cout << "allele:haplotype merit distribution" << std::endl;
//    for(auto& item: merit_dist){
//        std::cout << item.first << " : " << item.second << std::endl;
//    }

    //std::cout << "count of observed alleles: " << hostPool.getAlleleCounts() << std::endl;
//    std::vector<std::unordered_map<int, int>> allele_dists = hostPool.getAlleleDistributions();
//    std::cout << "Host allele distribution" << std::endl;
//    for(auto& dist : allele_dists){
//        for(auto& item : dist){
//            std::cout << item.first << " : " << item.second << std::endl;
//        }
//    }




    // i host generations with j pathogen generations each
    // each host totalHostGenerations contains
    //  n pathogen generations
    //      infection
    //      reproduction
    //      mutation
    //  reproduction
    //  mutation

    // begin pathogen generations per host generation
    for(int pathogen_generation = 0; pathogen_generation < config["infection"]["infections_per_generation"]; pathogen_generation++){
        totalPathogenGenerations++;

        // begin infection (occurs per pathogen generation)
        for(int host_species_index = 0; host_species_index < hostPool.hosts.size(); host_species_index++){
            unsigned long host_pop_size = hostPool.hosts[host_species_index].size();
            for(int host_index = 0; host_index < host_pop_size; host_index++) {
                Host& selectedHost = hostPool.hosts[host_species_index][host_index];

                for(int patho_species_index = 0; patho_species_index < pathogenPool.pathogens.size(); patho_species_index++){
                    int selectedPathogenIndex = rng.sampleIntUniUnsignedInt(0,pathogenPool.pathogens[patho_species_index].size() - 1);
                    Pathogen& selectedPathogen =  pathogenPool.pathogens[patho_species_index][selectedPathogenIndex];

                    // determine smallest lev distance of all host alleles to the selected pathogen haplotype
                    int smallest_lev = 99999;
                    for(const int& allele_id : selectedHost.chromosome_1_allele_ids){
                        const int lev_dist = meritCache.get(host_species_index, allele_id, patho_species_index, selectedPathogen.haplotype_id);
                        //std::cout << "lev dist: " << lev_dist << std::endl;
                        if(lev_dist < smallest_lev){
                            smallest_lev = lev_dist;
                        }
                    }
                    for(const int& allele_id : selectedHost.chromosome_2_allele_ids){
                        const int lev_dist = meritCache.get(host_species_index, allele_id, patho_species_index, selectedPathogen.haplotype_id);
                        //std::cout << "lev dist: " << lev_dist << std::endl;
                        if(lev_dist < smallest_lev){
                            smallest_lev = lev_dist;
                        }
                    }
                    //std::cout << "smallest found lev: " << smallest_lev << "\n";

                    // lev is lte than the provided threshold, at least one peptide of the pathogens haplotype is successfully presented by at least one mhc allele of the given host on either chromosome
                    if(smallest_lev <= config["infection"]["merit_threshold"]){
                        //std::cout << "successful presentation!" << std::endl;
                        selectedHost.antigen_presentation_count++;
                        selectedPathogen.no_infection_count++;
                    }else{
                        selectedHost.no_antigen_presentation_count++;
                        selectedPathogen.infection_count++;
                    }

                }
            }
        }

        // calculate fitness of hosts/pathogens and their sums after infection took place
        pathogenPool.updateFitness();

        // pathogen reproduction
        // infection count/fitness is cleared during this step!!
        double dice;
        int totalTries = 0;
        //TODO(JAN): find fast implementation for roulette wheel selection
        for(int patho_species_index = 0; patho_species_index < pathogenPool.pathogens.size(); patho_species_index++){
            int selectedPathogens = 0;

            // create vector of pathogens for this species
            std::vector<Pathogen> nextGenerationPathogens;
            nextGenerationPathogens.reserve(pathogenPool.pathogens[patho_species_index].size());

            unsigned int pathogenIdBase = pathogenPool.pathogens[patho_species_index].size() * totalPathogenGenerations;

            std::cout << "species total fitness: " << pathogenPool.fitness_sum[patho_species_index] << std::endl;
            while(selectedPathogens < pathogenPool.pathogens[patho_species_index].size()){
                dice = rng.sampleRealUniDouble(0, pathogenPool.fitness_sum[patho_species_index]);

                for(auto & pathogen : pathogenPool.pathogens[patho_species_index]){
                    totalTries++;
                    dice = dice - pathogen.fitness;
                    if(dice <= 0){
                        nextGenerationPathogens.emplace_back(Pathogen(pathogen.id, pathogenIdBase + selectedPathogens, config["pathogens"]["initial_fitness"], patho_species_index, pathogen.haplotype_id));
                        selectedPathogens++;
                        break;
                    }
                }
            }

            pathogenPool.pathogens[patho_species_index] = nextGenerationPathogens;
        }

        std::cout << "total tries: " << totalTries << std::endl;


        // get the allele distribution to skip cache filling for alleles that are not present in the population anymore (and therefore can't at any point in the future, too)
        std::vector<std::unordered_map<int, int>> host_allele_dist = hostPool.getAlleleDistributions();


        // pathogen mutation
        int haplotype_seq_length = config["pathogens"]["haplotype_sequence_length"];
        double mutation_rate_per_site = config["pathogens"]["mutation_rate_per_peptide"];
        std::string AS = config["aminoacids"];
        for(int patho_species_index = 0; patho_species_index < pathogenPool.pathogens.size(); patho_species_index++){
            for(Pathogen& currentPathogen : pathogenPool.pathogens[patho_species_index]){

                int mutationCount = rng.sampleBinomial(haplotype_seq_length, mutation_rate_per_site);
                if(mutationCount == 0) continue;

                std::cout << "mutation in " << currentPathogen.id << "\nold haplotype id: " << currentPathogen.haplotype_id << "\nsequence: " << pathogenAllelePool.alleles[patho_species_index][currentPathogen.haplotype_id].sequence << "\n";


                std::string newSequence = pathogenAllelePool.alleles[patho_species_index][currentPathogen.haplotype_id].sequence;

                for(int mutation_i = 0; mutation_i < mutationCount; mutation_i++){
                    unsigned int position = rng.sampleIntUniUnsignedInt(0, newSequence.size() - 1);
                    char newChar = AS[rng.sampleIntUniUnsignedInt(0, AS.size() -1)];
                    newSequence[position] = newChar;
                }

                unsigned long newHaplotypeId = pathogenAllelePool.addAllele(patho_species_index, newSequence);
                currentPathogen.haplotype_id = (int)newHaplotypeId;


                std::cout << "new haplotype id: " << newHaplotypeId << "\nnew sequence: " << newSequence << "\n";

                for( int host_species_i = 0; host_species_i < hostAllelePool.alleles.size(); host_species_i++ ){
                    for( auto &hostAllele : hostAllelePool.alleles[host_species_i] ){
                        //TODO(JAN): test this
                        if(host_allele_dist[host_species_i].find(hostAllele.id) != host_allele_dist[host_species_i].end()) continue; // skip this allele if no host actually carries it

                        int levDistance = Helper::generate_merit(hostAllele.sequence, newSequence);
                        meritCache.set(host_species_i, hostAllele.id, patho_species_index, newHaplotypeId, levDistance);
                    }
                }
            }
        }
    }

    // after n pathogen generations have passed (including n infections)
    // hosts move on in their lifecycle with reproduction and mutation

    // host reproduction
    //TODO(JAN): test all this
    std::cout << "Host reproduction\n";

    // infections are done, update and tally the fitness of the hosts
    hostPool.updateFitness();

    double dice;
    int totalTries = 0;

    for(int host_species_index = 0; host_species_index < hostPool.hosts.size(); host_species_index++){
        int selectedHosts = 0;
        int selectedParents = 0;

        // create vector of hosts for this species
        std::vector<Host> nextGenerationHosts;
        nextGenerationHosts.reserve(hostPool.hosts[host_species_index].size());

        unsigned int hostIdBase = hostPool.hosts[host_species_index].size() * totalHostGenerations;

        std::cout << "species total fitness: " << hostPool.fitness_sum[host_species_index] << std::endl;
        // need to select double the hosts -> two parents per new host
        while(selectedParents < hostPool.hosts[host_species_index].size() * 2){

            dice = rng.sampleRealUniDouble(0, hostPool.fitness_sum[host_species_index]);

            // select first parent
            int parentIndexFirst;
            for(int host_i = 0; host_i < hostPool.hosts[host_species_index].size(); host_i++){
                totalTries++;
                dice = dice - hostPool.hosts[host_species_index][host_i].fitness;
                if(dice <= 0){
                    parentIndexFirst = host_i;
                    selectedParents++;
                    break;
                }
            }

            // select second parent
            int parentIndexSecond;
            for(int host_i = 0; host_i < hostPool.hosts[host_species_index].size(); host_i++){
                totalTries++;
                dice = dice - hostPool.hosts[host_species_index][host_i].fitness;
                if(dice <= 0){
                    parentIndexSecond = host_i;
                    selectedParents++;
                    break;
                }
            }

            selectedHosts++;

            Host& parentFirst = hostPool.hosts[host_species_index][parentIndexFirst];
            Host& parentSecond = hostPool.hosts[host_species_index][parentIndexSecond];
            nextGenerationHosts.emplace_back(Host(parentFirst.id, parentSecond.id, hostIdBase + selectedHosts, config["hosts"]["initial_fitness"], host_species_index));
            Host& nextGenerationHost = nextGenerationHosts.back();

            if(rng.sampleRealUniFloat(0,1) < 0.5){
                nextGenerationHost.chromosome_1_allele_ids = parentFirst.chromosome_1_allele_ids;
            }else{
                nextGenerationHost.chromosome_1_allele_ids = parentFirst.chromosome_2_allele_ids;
            }

            if(rng.sampleRealUniFloat(0,1) < 0.5){
                nextGenerationHost.chromosome_2_allele_ids = parentSecond.chromosome_1_allele_ids;
            }else{
                nextGenerationHost.chromosome_2_allele_ids = parentSecond.chromosome_2_allele_ids;
            }

        }

        hostPool.hosts[host_species_index] = nextGenerationHosts;
    }


    // host mutation

    // get the allele distribution to skip cache filling for alleles that are not present in the population anymore (and therefore can't at any point in the future, too)
    std::vector<std::unordered_map<int, int>> pathogen_haplotype_dist = pathogenPool.getHaplotypeDistributions();

    int allele_seq_length = config["hosts"]["allele_sequence_length"];
    double mutation_rate_per_site = config["hosts"]["mutation_rate_per_peptide"];
    std::string AS = config["aminoacids"];
    for(int host_species_index = 0; host_species_index < hostPool.hosts.size(); host_species_index++){
        for(Host& currentHost : hostPool.hosts[host_species_index]){

            for(int& alleleId : currentHost.chromosome_1_allele_ids){
                int mutationCount = rng.sampleBinomial(allele_seq_length, mutation_rate_per_site);
                if(mutationCount == 0) continue;

                std::string newSequence = hostAllelePool.alleles[host_species_index][alleleId].sequence;

                std::cout << "mutation in allele: " << alleleId << "\nsequence: " << hostAllelePool.alleles[host_species_index][alleleId].sequence;

                for(int mutation_i = 0; mutation_i < mutationCount; mutation_i++){
                    unsigned int position = rng.sampleIntUniUnsignedInt(0, newSequence.size() - 1);
                    char newChar = AS[rng.sampleIntUniUnsignedInt(0, AS.size() -1)];
                    newSequence[position] = newChar;
                }

                unsigned long newAlleleId = hostAllelePool.addAllele(host_species_index, newSequence);
                alleleId = (int)newAlleleId;


                std::cout << "new allele id: " << newAlleleId << "\nnew sequence: " << newSequence << "\n";

                for( int pathogen_species_i = 0; pathogen_species_i < pathogenAllelePool.alleles.size(); pathogen_species_i++ ){
                    for( auto &pathogenHaplotype : pathogenAllelePool.alleles[pathogen_species_i] ){
                        //TODO(JAN): test this
                        if(pathogen_haplotype_dist[pathogen_species_i].find(pathogenHaplotype.id) != pathogen_haplotype_dist[pathogen_species_i].end()) continue; // skip this allele if no host actually carries it

                        int levDistance = Helper::generate_merit(pathogenHaplotype.sequence, newSequence);
                        meritCache.set(host_species_index, newAlleleId, pathogen_species_i, pathogenHaplotype.id, levDistance);
                    }
                }

            }

            for(int& alleleId : currentHost.chromosome_2_allele_ids){
                int mutationCount = rng.sampleBinomial(allele_seq_length, mutation_rate_per_site);
                if(mutationCount == 0) continue;

                std::string newSequence = hostAllelePool.alleles[host_species_index][alleleId].sequence;

                std::cout << "mutation in allele: " << alleleId << "\nsequence: " << hostAllelePool.alleles[host_species_index][alleleId].sequence;

                for(int mutation_i = 0; mutation_i < mutationCount; mutation_i++){
                    unsigned int position = rng.sampleIntUniUnsignedInt(0, newSequence.size() - 1);
                    char newChar = AS[rng.sampleIntUniUnsignedInt(0, AS.size() -1)];
                    newSequence[position] = newChar;
                }

                unsigned long newAlleleId = hostAllelePool.addAllele(host_species_index, newSequence);
                alleleId = (int)newAlleleId;


                std::cout << "new allele id: " << newAlleleId << "\nnew sequence: " << newSequence << "\n";

                for( int pathogen_species_i = 0; pathogen_species_i < pathogenAllelePool.alleles.size(); pathogen_species_i++ ){
                    for( auto &pathogenHaplotype : pathogenAllelePool.alleles[pathogen_species_i] ){
                        //TODO(JAN): test this
                        if(pathogen_haplotype_dist[pathogen_species_i].find(pathogenHaplotype.id) != pathogen_haplotype_dist[pathogen_species_i].end()) continue; // skip this allele if no host actually carries it

                        int levDistance = Helper::generate_merit(pathogenHaplotype.sequence, newSequence);
                        meritCache.set(host_species_index, newAlleleId, pathogen_species_i, pathogenHaplotype.id, levDistance);
                    }
                }

            }

        }
    }



    std::cout << "generation " << totalHostGenerations << " complete\n";

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