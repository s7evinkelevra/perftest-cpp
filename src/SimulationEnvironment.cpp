//
// Created by Jan on 20.04.2022.
//

#include "SimulationEnvironment.h"

#include <iostream>
#include <utility>
#include <filesystem>
#include <thread>
#include <omp.h>


SimulationEnvironment::SimulationEnvironment(json initialConfig){
    config = std::move(initialConfig);
    totalHostGenerations = 0;
    totalPathogenGenerations = 0;

    bInfection = true;
    bHostFitnessproportionalReproduction = true;
    bPathogenFitnessproportionalReproduction = true;
    bHostMutation = true;
    bPathogenMutation = true;

    rng = Random();

    unsigned int threads_wanted = config["thread_n"];
    unsigned int threads_available = std::thread::hardware_concurrency();
    thread_count = std::min({threads_available, threads_wanted});
    std::cout << "threads wanted: " << threads_wanted << "\nthreads available: " << threads_available << "\n";
    omp_set_num_threads(thread_count);
}

void SimulationEnvironment::initializeHostAllelePool() {
    hostAllelePool.alleles.resize(config["hosts"]["species_n"]);
    hostAllelePool.total_allele_counts_per_species.resize(config["hosts"]["species_n"]);
    for( int species_i = 0; species_i < config["hosts"]["species_n"]; species_i ++ ){

        hostAllelePool.alleles[species_i].reserve(config["hosts"]["alleles_per_species_initial"]);
        hostAllelePool.total_allele_counts_per_species[species_i] = config["hosts"]["alleles_per_species_initial"];
        for(int i = 0; i < config["hosts"]["alleles_per_species_initial"]; i++){
            hostAllelePool.alleles[species_i].emplace(i, Allele(-1, i, totalHostGenerations, generateSequence(config["hosts"]["allele_sequence_length"])));
        }
    }
}

void SimulationEnvironment::initializePathogenAllelePool() {
    pathogenAllelePool.alleles.resize(config["pathogens"]["species_n"]);
    pathogenAllelePool.total_allele_counts_per_species.resize(config["pathogens"]["species_n"]);
    for( int species_i = 0; species_i < config["pathogens"]["species_n"]; species_i++ ){

        pathogenAllelePool.alleles[species_i].reserve(config["pathogens"]["haplotypes_per_species_initial"]);
        pathogenAllelePool.total_allele_counts_per_species[species_i] = config["pathogens"]["haplotypes_per_species_initial"];
        for(int i = 0; i < config["pathogens"]["haplotypes_per_species_initial"]; i++){
            pathogenAllelePool.alleles[species_i].emplace(i, Allele(-1, i, totalPathogenGenerations, generateSequence(config["pathogens"]["haplotype_sequence_length"])));
        }
    }
}

void SimulationEnvironment::initializeMeritCache() {
    for( int host_species_i = 0; host_species_i < config["hosts"]["species_n"]; host_species_i++ ){
        for( auto & [host_allele_id, hostAllele] : hostAllelePool.alleles[host_species_i] ){
            for( int patho_species_i = 0; patho_species_i < config["pathogens"]["species_n"]; patho_species_i++ ){
                for( auto &[pathogen_allele_id, pathogenAllele] : pathogenAllelePool.alleles[patho_species_i] ) {
                    int levDistance = Helper::generate_merit(hostAllele.sequence, pathogenAllele.sequence);
                    meritCache.set(host_species_i, hostAllele.id, patho_species_i, pathogenAllele.id, levDistance);
                }
            }
        }
    }
};

void SimulationEnvironment::initializeHostPool() {
    hostPool.max_loci_count = 1;
    hostPool.hosts.resize(config["hosts"]["species_n"]);
    hostPool.fitness_sum.resize(config["hosts"]["species_n"]);
    int initialFitness = config["hosts"]["initial_fitness"];

    for( int species_i = 0; species_i < config["hosts"]["species_n"]; species_i ++) {
        hostPool.hosts[species_i].reserve(config["hosts"]["n"]);
        for( int i = 0; i < config["hosts"]["n"]; i++ ){
            hostPool.hosts[species_i].emplace_back(Host(0, initialFitness, 0, initialFitness, i, initialFitness, species_i));
            for(int j = 0; j < config["hosts"]["genes_per_chromosome_initial"]; j++){
                auto& allelePool = hostAllelePool.alleles[species_i];

//                auto randomAlleleId_1_it = std::next(allelePool.begin(), rng.sampleIntUniUnsignedInt(0, allelePool.size() - 1));
//                auto randomAlleleId_2_it = std::next(allelePool.begin(), rng.sampleIntUniUnsignedInt(0, allelePool.size() - 1));

                int randomAlleleId_1 = rng.sampleIntUniUnsignedInt(0, allelePool.size() - 1);
                int randomAlleleId_2 = rng.sampleIntUniUnsignedInt(0, allelePool.size() - 1);

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
            auto& allelePool = pathogenAllelePool.alleles[species_i];

            //auto randomHaplotypeId_it = std::next(allelePool.begin(), rng.sampleIntUniUnsignedInt(0, pathogenAllelePool.alleles[species_i].size() - 1));
            int randomHaplotypeId = rng.sampleIntUniUnsignedInt(0, pathogenAllelePool.alleles[species_i].size() - 1);


            pathogenPool.pathogens[species_i].emplace_back(Pathogen(0,i,initialFitness,species_i,randomHaplotypeId));
        }
    }
}

void SimulationEnvironment::initializeOutputFiles() {
    // create CSV writers and files
    std::string output_base_path = config["output"]["output_path"];

    if(std::filesystem::is_directory(output_base_path) || std::filesystem::exists(output_base_path)){
        std::filesystem::remove_all(output_base_path);
    }

    std::filesystem::create_directory(output_base_path);

    // copy config file to output dir
    std::ofstream configFileStream(output_base_path + "config.json");
    configFileStream << config.dump(4);
    configFileStream.close();

    std::vector allele_CSV_headers = {"generation", "species", "locus_id", "parent_id", "allele_id","created_at", "count", "frequency"};
    std::vector locus_CSV_headers = {"generation", "species", "locus_id", "allele_count", "allelic_richness", "H_e", "H_o", "HWE"};


    std::vector host_CSV_headers = {"generation", "species", "id", "parent_1_id", "parent_2_id", "successful_presentations", "unsuccessful_presentations", "total_presentations", "fitness"};
    hostDataCSV = std::make_unique<CSVWriter>(output_base_path + "host_data.csv", ";");
    hostDataCSV->addRow(host_CSV_headers.begin(), host_CSV_headers.end());

    std::vector host_genome_CSV_headers = {"generation", "species", "id", "locus_id", "allele_1_id", "allele_2_id"};
    hostGenomeDataCSV = std::make_unique<CSVWriter>(output_base_path + "host_genome_data.csv", ";");
    hostGenomeDataCSV->addRow(host_genome_CSV_headers.begin(), host_genome_CSV_headers.end());

    hostAlleleDataCSV = std::make_unique<CSVWriter>(output_base_path + "host_allele_data.csv", ";");
    hostAlleleDataCSV->addRow(allele_CSV_headers.begin(), allele_CSV_headers.end());

    hostLocusDataCSV = std::make_unique<CSVWriter>(output_base_path + "host_locus_data.csv", ";");
    hostLocusDataCSV->addRow(locus_CSV_headers.begin(), locus_CSV_headers.end());


    std::vector pathogen_CSV_headers = {"generation", "species", "id", "parent_id", "successful_infections", "unsuccessful_infections", "total_infections", "fitness"};
    pathogenDataCSV = std::make_unique<CSVWriter>(output_base_path + "pathogen_data.csv", ";");
    pathogenDataCSV->addRow(pathogen_CSV_headers.begin(), pathogen_CSV_headers.end());

    std::vector pathogen_genome_CSV_headers = {"generation", "species", "id", "haplotype_id"};
    pathogenGenomeDataCSV = std::make_unique<CSVWriter>(output_base_path + "pathogen_genome_data.csv", ";");
    pathogenGenomeDataCSV->addRow(pathogen_genome_CSV_headers.begin(), pathogen_genome_CSV_headers.end());

    pathogenAlleleDataCSV = std::make_unique<CSVWriter>(output_base_path + "pathogen_allele_data.csv", ";");
    pathogenAlleleDataCSV->addRow(allele_CSV_headers.begin(), allele_CSV_headers.end());

    pathogenLocusDataCSV = std::make_unique<CSVWriter>(output_base_path + "pathogen_locus_data.csv", ";");
    pathogenLocusDataCSV->addRow(locus_CSV_headers.begin(), locus_CSV_headers.end());

    std::vector meta_CSV_headers = {"generation", "pathogen_generation", "bInfection", "bHostFitnessproportionalReproduction", "bPathogenFitnessproportionalReproduction", "bHostMutation", "bPathogenMutation", "step_duration"};
    metaDataCSV = std::make_unique<CSVWriter>(output_base_path + "meta_data.csv", ";");
    metaDataCSV->addRow(meta_CSV_headers.begin(), meta_CSV_headers.end());

}


void SimulationEnvironment::printHost(int species, int index){
    Host& host = hostPool.hosts[species][index];
    host.print();

    std::cout << "  chromosome 1 alleles: " << std::endl;
    for(const auto &alleleId : host.chromosome_1_allele_ids){
        Allele& allele = hostAllelePool.alleles[species].at(alleleId);
        std::cout << "  allele id: " << allele.id << " -> " << allele.sequence << std::endl;
    }

    std::cout << "  chromosome 2 alleles: " << std::endl;
    for(const auto &alleleId : host.chromosome_2_allele_ids){
        Allele& allele = hostAllelePool.alleles[species].at(alleleId);
        std::cout << "  allele id: " << allele.id << " -> " << allele.sequence << std::endl;
    }
}

void SimulationEnvironment::printPathogen(int species, int index){
    Pathogen& pathogen = pathogenPool.pathogens[species][index];
    pathogen.print();

    Allele& haplotype = pathogenAllelePool.alleles[species].at(pathogen.haplotype_id);
    std::cout << "  haplotype id: " << haplotype.id << " -> " << haplotype.sequence << std::endl;
}


void SimulationEnvironment::initialize() {
    initializeHostAllelePool();
    initializePathogenAllelePool();

    // after burn-in, none of the originial alleles are left anyways...
    initializeMeritCache();

    initializeHostPool();
    initializePathogenPool();

    initializeOutputFiles();
    writeAllData();
}

// implement single simulation step
// this represents a single host totalHostGenerations with n pathogen generations
void SimulationEnvironment::step() {
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

    lastStepStart = std::chrono::steady_clock::now();

    hostGeneration();
    lastStepEnd = std::chrono::steady_clock::now();


    int allele_data_interval = config["output"]["allele_data_interval"];
    if(totalHostGenerations % allele_data_interval == 0){
        writeHostAlleleData();
        writePathogenAlleleData();
    }

    int locus_data_interval = config["output"]["locus_data_interval"];
    if(totalHostGenerations % locus_data_interval == 0){
        writeHostLocusData();
        writePathogenLocusData();
    }

    int meta_data_interval = config["output"]["meta_data_interval"];
    if(totalHostGenerations % meta_data_interval == 0){
        writeMetaData();
    }



    std::cout << "generation " << totalHostGenerations << " complete\n";
    std::cout << "time elapsed: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(lastStepEnd - lastStepStart).count()
              << " ms" << std::endl;
}

void SimulationEnvironment::hostGeneration() {
    totalHostGenerations++;
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
        pathogenGeneration();
    }

    // after n pathogen generations have passed (including n infections)
    // hosts move on in their lifecycle with reproduction and mutation

    // infections are done, update and tally the fitness of the hosts
    // necessary for logging
    hostPool.updateFitness();

    int individual_data_interval = config["output"]["individual_data_interval"];
    if(totalHostGenerations % individual_data_interval == 0){
        writeHostData();
        writeHostGenomeData();
    }

    // host reproduction
    //TODO(JAN): test all this
    std::cout << "Host reproduction\n";
    if(bHostFitnessproportionalReproduction){
        hostReproduction();
    }else{
        hostReproductionRandom();
    }

    // host mutation
    if(bHostMutation){
        hostMutation();
    }
}

void SimulationEnvironment::hostMutation() {
    // get the allele distribution to skip cache filling for alleles that are not present in the population anymore (and therefore can't at any point in the future, too)
    //std::vector<std::unordered_map<int, int>> pathogen_haplotype_dist = pathogenPool.getHaplotypeDistributions();

    int allele_seq_length = config["hosts"]["allele_sequence_length"];
    double mutation_rate_per_site = config["hosts"]["mutation_rate_per_peptide"];
    std::string AS = config["aminoacids"];
    for(int host_species_index = 0; host_species_index < hostPool.hosts.size(); host_species_index++){
        for(Host& currentHost : hostPool.hosts[host_species_index]){

            for(int& alleleId : currentHost.chromosome_1_allele_ids){
                int mutationCount = rng.sampleBinomial(allele_seq_length, mutation_rate_per_site);
                if(mutationCount == 0) continue;

                std::string newSequence = hostAllelePool.alleles[host_species_index].at(alleleId).sequence;

                //std::cout << "mutation in allele: " << alleleId << "\nsequence: " << hostAllelePool.alleles[host_species_index][alleleId].sequence << "\n";

                for(int mutation_i = 0; mutation_i < mutationCount; mutation_i++){
                    unsigned int position = rng.sampleIntUniUnsignedInt(0, newSequence.size() - 1);
                    char newChar = AS[rng.sampleIntUniUnsignedInt(0, AS.size() - 1)];
                    newSequence[position] = newChar;
                }

                unsigned long newAlleleId = hostAllelePool.addAllele(host_species_index, alleleId, totalHostGenerations, newSequence);
                alleleId = (int)newAlleleId;
            }

            for(int& alleleId : currentHost.chromosome_2_allele_ids){
                int mutationCount = rng.sampleBinomial(allele_seq_length, mutation_rate_per_site);
                if(mutationCount == 0) continue;

                std::string newSequence = hostAllelePool.alleles[host_species_index].at(alleleId).sequence;

                //std::cout << "mutation in allele: " << alleleId << "\nsequence: " << hostAllelePool.alleles[host_species_index][alleleId].sequence << "\n";

                for(int mutation_i = 0; mutation_i < mutationCount; mutation_i++){
                    unsigned int position = rng.sampleIntUniUnsignedInt(0, newSequence.size() - 1);
                    char newChar = AS[rng.sampleIntUniUnsignedInt(0, AS.size() - 1)];
                    newSequence[position] = newChar;
                }

                unsigned long newAlleleId = hostAllelePool.addAllele(host_species_index, alleleId, totalHostGenerations, newSequence);
                alleleId = (int)newAlleleId;

            }

        }
    }
}

void SimulationEnvironment::hostReproduction() {
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

            // select first parent
            dice = rng.sampleRealUniDouble(0, hostPool.fitness_sum[host_species_index]);
            int parent_1_id;
            for(int host_i = 0; host_i < hostPool.hosts[host_species_index].size(); host_i++){
                totalTries++;
                dice = dice - hostPool.hosts[host_species_index][host_i].fitness;
                if(dice <= 0){
                    parent_1_id = host_i;
                    selectedParents++;
                    break;
                }
            }

            // select second parent
            dice = rng.sampleRealUniDouble(0, hostPool.fitness_sum[host_species_index]);
            int parent_2_id;
            for(int host_i = 0; host_i < hostPool.hosts[host_species_index].size(); host_i++){
                totalTries++;
                dice = dice - hostPool.hosts[host_species_index][host_i].fitness;
                if(dice <= 0){
                    parent_2_id = host_i;
                    selectedParents++;
                    break;
                }
            }

            selectedHosts++;

            Host& parent_1 = hostPool.hosts[host_species_index][parent_1_id];
            Host& parent_2 = hostPool.hosts[host_species_index][parent_2_id];
            nextGenerationHosts.emplace_back(Host(parent_1.id, parent_1.fitness, parent_2.id, parent_2.fitness, hostIdBase + selectedHosts, config["hosts"]["initial_fitness"], host_species_index));
            Host& nextGenerationHost = nextGenerationHosts.back();

            if(rng.sampleRealUniFloat(0, 1) < 0.5){
                nextGenerationHost.chromosome_1_allele_ids = parent_1.chromosome_1_allele_ids;
            }else{
                nextGenerationHost.chromosome_1_allele_ids = parent_1.chromosome_2_allele_ids;
            }

            if(rng.sampleRealUniFloat(0, 1) < 0.5){
                nextGenerationHost.chromosome_2_allele_ids = parent_2.chromosome_1_allele_ids;
            }else{
                nextGenerationHost.chromosome_2_allele_ids = parent_2.chromosome_2_allele_ids;
            }

        }

        hostPool.hosts[host_species_index] = nextGenerationHosts;
    }
}

//TODO(JAN): test this
void SimulationEnvironment::hostReproductionRandom() {
    for(int host_species_index = 0; host_species_index < hostPool.hosts.size(); host_species_index++){

        unsigned int host_pop_size = hostPool.hosts[host_species_index].size();
        // create vector of hosts for this species
        std::vector<Host> nextGenerationHosts;
        nextGenerationHosts.reserve(host_pop_size);

        unsigned int hostIdBase = hostPool.hosts[host_species_index].size() * totalHostGenerations;

        for(int new_host_i = 0; new_host_i < host_pop_size; new_host_i++){
            unsigned int parent_1_id = rng.sampleIntUniUnsignedInt(0, host_pop_size - 1);
            unsigned int parent_2_id = rng.sampleIntUniUnsignedInt(0, host_pop_size - 1);

            Host& parent_1 = hostPool.hosts[host_species_index][parent_1_id];
            Host& parent_2 = hostPool.hosts[host_species_index][parent_2_id];

            nextGenerationHosts.emplace_back(Host(parent_1_id, parent_1.fitness, parent_2_id, parent_2.fitness, hostIdBase + new_host_i, config["hosts"]["initial_fitness"], host_species_index));
            Host& nextGenerationHost = nextGenerationHosts.back();


            if(rng.sampleRealUniFloat(0,1) < 0.5){
                nextGenerationHost.chromosome_1_allele_ids = parent_1.chromosome_1_allele_ids;
            }else{
                nextGenerationHost.chromosome_1_allele_ids = parent_1.chromosome_2_allele_ids;
            }

            if(rng.sampleRealUniFloat(0,1) < 0.5){
                nextGenerationHost.chromosome_2_allele_ids = parent_2.chromosome_1_allele_ids;
            }else{
                nextGenerationHost.chromosome_2_allele_ids = parent_2.chromosome_2_allele_ids;
            }
        }
        hostPool.hosts[host_species_index] = nextGenerationHosts;
    }
}


void SimulationEnvironment::pathogenGeneration() {
    totalPathogenGenerations++;

    // begin infection (occurs per pathogen generation)
    if(bInfection){
        infection();
    }

    // calculate fitness of hosts/pathogens and their sums after infection took place
    pathogenPool.updateFitness();

    int individual_data_interval = config["output"]["individual_data_interval"];
    int pathogen_generations_per_host_generation = config["infection"]["infections_per_generation"];
    if((totalHostGenerations % (individual_data_interval * pathogen_generations_per_host_generation)) == 0){
        writePathogenData();
        writePathogenGenomeData();
    }

    if(bPathogenFitnessproportionalReproduction){
        pathogenReproduction();
    }else{
        pathogenReproductionRandom();
    }

    if(bPathogenMutation){
        pathogenMutation();
    }
}

void SimulationEnvironment::pathogenMutation() {
    // get the allele distribution to skip cache filling for alleles that are not present in the population anymore (and therefore can't at any point in the future, too)
    //std::vector<std::unordered_map<int, int>> host_allele_dist_across_all_loci = hostPool.getAlleleDistributionAcrossAllLociPerSpecies();


    // pathogen mutation
    int haplotype_seq_length = config["pathogens"]["haplotype_sequence_length"];
    double mutation_rate_per_site = config["pathogens"]["mutation_rate_per_peptide"];
    std::string AS = config["aminoacids"];
    for(int patho_species_index = 0; patho_species_index < pathogenPool.pathogens.size(); patho_species_index++){
        for(Pathogen& currentPathogen : pathogenPool.pathogens[patho_species_index]){

            int mutationCount = rng.sampleBinomial(haplotype_seq_length, mutation_rate_per_site);
            if(mutationCount == 0) continue;

            //std::cout << "mutation in " << currentPathogen.id << "\nold haplotype id: " << currentPathogen.haplotype_id << "\nsequence: " << pathogenAllelePool.alleles[patho_species_index][currentPathogen.haplotype_id].sequence << "\n";


            std::string newSequence = pathogenAllelePool.alleles[patho_species_index].at(currentPathogen.haplotype_id).sequence;

            for(int mutation_i = 0; mutation_i < mutationCount; mutation_i++){
                unsigned int position = rng.sampleIntUniUnsignedInt(0, newSequence.size() - 1);
                char newChar = AS[rng.sampleIntUniUnsignedInt(0, AS.size() - 1)];
                newSequence[position] = newChar;
            }

            unsigned long newHaplotypeId = pathogenAllelePool.addAllele(patho_species_index, currentPathogen.haplotype_id, totalPathogenGenerations, newSequence);
            currentPathogen.haplotype_id = (int)newHaplotypeId;


            //std::cout << "new haplotype id: " << newHaplotypeId << "\nnew sequence: " << newSequence << "\n";
            /*
            for(int host_species_i = 0; host_species_i < hostAllelePool.alleles.size(); host_species_i++ ){
                for( auto &hostAllele : hostAllelePool.alleles[host_species_i] ){
                    //TODO(JAN): test this
                    if(host_allele_dist_across_all_loci[host_species_i].find(hostAllele.id) != host_allele_dist_across_all_loci[host_species_i].end()) continue; // skip this allele if no host actually carries it

                    int levDistance = Helper::generate_merit(hostAllele.sequence, newSequence);
                    meritCache.set(host_species_i, hostAllele.id, patho_species_index, (int)newHaplotypeId, levDistance);
                }
            }
            */
        }
    }
}

void SimulationEnvironment::pathogenReproduction() {


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
}

void SimulationEnvironment::pathogenReproductionRandom() {
    for(int patho_species_index = 0; patho_species_index < pathogenPool.pathogens.size(); patho_species_index++){

        // create vector of pathogens for this species
        std::vector<Pathogen> nextGenerationPathogens;
        unsigned int pathogen_pop_size = pathogenPool.pathogens[patho_species_index].size();
        nextGenerationPathogens.reserve(pathogen_pop_size);

        unsigned int pathogenIdBase = pathogenPool.pathogens[patho_species_index].size() * totalPathogenGenerations;


        for(int new_pathogen_i = 0; new_pathogen_i < pathogen_pop_size; new_pathogen_i++){
            int parent_id = rng.sampleIntUniUnsignedInt(0, pathogen_pop_size - 1);
            Pathogen& pathogen_parent = pathogenPool.pathogens[patho_species_index][parent_id];
            nextGenerationPathogens.emplace_back(Pathogen(parent_id, pathogenIdBase + new_pathogen_i, config["pathogens"]["initial_fitness"], patho_species_index, pathogen_parent.haplotype_id));
        }

        pathogenPool.pathogens[patho_species_index] = nextGenerationPathogens;
    }


}

void SimulationEnvironment::infection() {
    for(int host_species_index = 0; host_species_index < hostPool.hosts.size(); host_species_index++){
        unsigned long host_pop_size = hostPool.hosts[host_species_index].size();
        for(int host_index = 0; host_index < host_pop_size; host_index++) {
            Host& selectedHost = hostPool.hosts[host_species_index][host_index];

            for(int patho_species_index = 0; patho_species_index < pathogenPool.pathogens.size(); patho_species_index++){
                int selectedPathogenIndex = rng.sampleIntUniUnsignedInt(0, pathogenPool.pathogens[patho_species_index].size() - 1);
                Pathogen& selectedPathogen =  pathogenPool.pathogens[patho_species_index][selectedPathogenIndex];

                // determine smallest lev distance of all host alleles to the selected pathogen haplotype
                int smallest_lev = 99999;
                for(const int& allele_id : selectedHost.chromosome_1_allele_ids){
                    //TODO(JAN): instead of generating merit values for each new allele:haplotype combination, check here if combination is present, else generate and set it
                    // -> fill cache on demand
                    int lev_dist = 9999;
                    if(meritCache.exists(host_species_index, allele_id, patho_species_index, selectedPathogen.haplotype_id)){
                        lev_dist = meritCache.get(host_species_index, allele_id, patho_species_index, selectedPathogen.haplotype_id);
                    }else{
                        Allele& hostAllele = hostAllelePool.alleles[host_species_index].at(allele_id);
                        Allele& pathogenHaplotype = pathogenAllelePool.alleles[patho_species_index].at(selectedPathogen.haplotype_id);
                        lev_dist = Helper::generate_merit(hostAllele.sequence, pathogenHaplotype.sequence);
                        meritCache.set(host_species_index, allele_id, patho_species_index, selectedPathogen.haplotype_id, lev_dist);
                    }
                    //std::cout << "lev dist: " << lev_dist << std::endl;
                    if(lev_dist < smallest_lev){
                        smallest_lev = lev_dist;
                    }
                }
                for(const int& allele_id : selectedHost.chromosome_2_allele_ids){
                    int lev_dist = 9999;
                    if(meritCache.exists(host_species_index, allele_id, patho_species_index, selectedPathogen.haplotype_id)){
                        lev_dist = meritCache.get(host_species_index, allele_id, patho_species_index, selectedPathogen.haplotype_id);
                    }else{
                        Allele& hostAllele = hostAllelePool.alleles[host_species_index].at(allele_id);
                        Allele& pathogenHaplotype = pathogenAllelePool.alleles[patho_species_index].at(selectedPathogen.haplotype_id);
                        lev_dist = Helper::generate_merit(hostAllele.sequence, pathogenHaplotype.sequence);
                        meritCache.set(host_species_index, allele_id, patho_species_index, selectedPathogen.haplotype_id, lev_dist);
                    }
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
}

std::string SimulationEnvironment::generateSequence(int length) {
    std::string AS = config["aminoacids"];
    std::string tmp_s;
    tmp_s.reserve(length);

    for (int i = 0; i < length; ++i) {
        tmp_s += AS[rng.sampleIntUniUnsignedInt(0, AS.length() - 1)];
    }

    return tmp_s;
}

void SimulationEnvironment::setDefaultMode() {
    bInfection = true;
    bHostFitnessproportionalReproduction = true;
    bPathogenFitnessproportionalReproduction = true;
    bHostMutation = true;
    bPathogenMutation = true;
}

void SimulationEnvironment::setBurnInMode() {
    bInfection = false;
    bHostFitnessproportionalReproduction = false;
    bPathogenFitnessproportionalReproduction = false;
    bHostMutation = true;
    bPathogenMutation = true;
}

void SimulationEnvironment::setNoCoevolutionMode() {
    bInfection = true;
    bHostFitnessproportionalReproduction = true;
    bPathogenFitnessproportionalReproduction = false;
    bHostMutation = true;
    bPathogenMutation = true;
}

void SimulationEnvironment::writeHostData() {
    std::cout << "writing host data" << "\n";

    for(int host_species_i = 0; host_species_i < hostPool.hosts.size(); host_species_i++){
        for(int host_i = 0; host_i < hostPool.hosts[host_species_i].size(); host_i++){

            Host& currentHost = hostPool.hosts[host_species_i][host_i];
            std::vector<std::string> props = {
                    std::to_string(totalHostGenerations),
                    std::to_string(host_species_i),
                    std::to_string(currentHost.id),
                    std::to_string(currentHost.parent_id_1),
                    std::to_string(currentHost.parent_id_2),
                    std::to_string(currentHost.antigen_presentation_count),
                    std::to_string(currentHost.no_antigen_presentation_count),
                    std::to_string(currentHost.antigen_presentation_count + currentHost.no_antigen_presentation_count),
                    std::to_string(currentHost.fitness)
                    };
            hostDataCSV->addRow(props.begin(), props.end());
        }
    }
}

void SimulationEnvironment::writeHostGenomeData() {
    std::cout << "writing host genome data" << "\n";

    for(int host_species_i = 0; host_species_i < hostPool.hosts.size(); host_species_i++){
        for(int host_i = 0; host_i < hostPool.hosts[host_species_i].size(); host_i++){
            Host& currentHost = hostPool.hosts[host_species_i][host_i];

            for(int locus_i = 0; locus_i < currentHost.chromosome_1_allele_ids.size(); locus_i++){
                std::vector<std::string> props = {
                        std::to_string(totalHostGenerations),
                        std::to_string(host_species_i),
                        std::to_string(currentHost.id),
                        std::to_string(locus_i),
                        std::to_string(currentHost.chromosome_1_allele_ids[locus_i]),
                        std::to_string(currentHost.chromosome_2_allele_ids[locus_i]),
                };
                hostGenomeDataCSV->addRow(props.begin(), props.end());
            }

        }
    }
}

void SimulationEnvironment::writeHostAlleleData() {
    std::cout << "writing host allele data" << "\n";
    for(int host_species_i = 0; host_species_i < hostPool.hosts.size(); host_species_i++){

        for(int locus_i = 0; locus_i < hostPool.max_loci_count; locus_i++){
            std::unordered_map<int,int> current_allele_dist = hostPool.getAlleleDistribution(host_species_i, locus_i);
            int allele_count_sum = 0;
            for(auto& item : current_allele_dist){
                allele_count_sum += item.second;
            }

            for(auto& item : current_allele_dist) {
                Allele &current_allele = hostAllelePool.alleles[host_species_i].at(item.first);

                std::vector<std::string> props = {
                        std::to_string(totalHostGenerations),
                        std::to_string(host_species_i),
                        std::to_string(locus_i),
                        std::to_string(current_allele.parentId),
                        std::to_string(item.first),
                        std::to_string(current_allele.createdAtGeneration),
                        std::to_string(item.second),
                        std::to_string((float) item.second / (float) allele_count_sum)
                };
                hostAlleleDataCSV->addRow(props.begin(), props.end());
            }
        }
    }

}

void SimulationEnvironment::writeMetaData() {
    std::cout << "writing meta data" << "\n";

    std::vector<std::string> props = {
            std::to_string(totalHostGenerations),
            std::to_string(totalPathogenGenerations),
            std::to_string(bInfection),
            std::to_string(bHostFitnessproportionalReproduction),
            std::to_string(bPathogenFitnessproportionalReproduction),
            std::to_string(bHostMutation),
            std::to_string(bPathogenMutation),
            std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(lastStepEnd - lastStepStart).count())
    };
    metaDataCSV->addRow(props.begin(), props.end());
}

void SimulationEnvironment::writePathogenData() {
    std::cout << "writing pathogen data" << "\n";

    for (int patho_species_i = 0; patho_species_i < pathogenPool.pathogens.size(); patho_species_i++){
        for (int pathogen_i = 0; pathogen_i < pathogenPool.pathogens[patho_species_i].size(); pathogen_i++){

            Pathogen &currentPathogen = pathogenPool.pathogens[patho_species_i][pathogen_i];
            std::vector<std::string> props = {
                std::to_string(totalPathogenGenerations),
                std::to_string(patho_species_i),
                std::to_string(currentPathogen.id),
                std::to_string(currentPathogen.parent_id),
                std::to_string(currentPathogen.infection_count),
                std::to_string(currentPathogen.no_infection_count),
                std::to_string(currentPathogen.infection_count + currentPathogen.no_infection_count),
                std::to_string(currentPathogen.fitness)};
            pathogenDataCSV->addRow(props.begin(), props.end());
        }
    }
}

void SimulationEnvironment::writePathogenGenomeData() {
    std::cout << "writing pathogen genome data" << "\n";

    for (int patho_species_i = 0; patho_species_i < pathogenPool.pathogens.size(); patho_species_i++){
        for (int pathogen_i = 0; pathogen_i < pathogenPool.pathogens[patho_species_i].size(); pathogen_i++){

            Pathogen &currentPathogen = pathogenPool.pathogens[patho_species_i][pathogen_i];
            std::vector<std::string> props = {
                std::to_string(totalPathogenGenerations),
                std::to_string(patho_species_i),
                std::to_string(currentPathogen.id),
                std::to_string(currentPathogen.haplotype_id)};
            pathogenGenomeDataCSV->addRow(props.begin(), props.end());
        }
    }
}

void SimulationEnvironment::writePathogenAlleleData() {
    std::cout << "writing pathogen allele data" << "\n";

    for(int patho_species_i = 0; patho_species_i < pathogenPool.pathogens.size(); patho_species_i++){
        std::unordered_map<int,int> haplotype_dist = pathogenPool.getHaplotypeDistribution(patho_species_i);
        int haplotype_count_sum = 0;
        for(auto& item : haplotype_dist){
            haplotype_count_sum += item.second;
        }


        for(auto& item : haplotype_dist){
            Allele& currentHaplotype = pathogenAllelePool.alleles[patho_species_i].at(item.first);

            std::vector<std::string> props = {
                std::to_string(totalPathogenGenerations),
                std::to_string(patho_species_i),
                std::to_string(0),
                std::to_string(currentHaplotype.parentId),
                std::to_string(item.first),
                std::to_string(currentHaplotype.createdAtGeneration),
                std::to_string(item.second),
                std::to_string((float)item.second/(float)haplotype_count_sum)
            };

            pathogenAlleleDataCSV->addRow(props.begin(), props.end());
        }
    }
}

void SimulationEnvironment::writeHostLocusData() {
    std::cout << "writing host locus data [not implemented]" << "\n";

}

void SimulationEnvironment::writePathogenLocusData() {
    std::cout << "writing pathogen locus data [not implemented]" << "\n";

}

void SimulationEnvironment::writeAllData() {
    writeHostData();
    writeHostGenomeData();

    writePathogenData();
    writePathogenGenomeData();

    writeHostAlleleData();
    writePathogenAlleleData();

    writeHostLocusData();
    writePathogenLocusData();

    writeMetaData();
}


void SimulationEnvironment::purgeUnusedAlleles(){
    for(int host_species_i = 0; host_species_i < hostAllelePool.alleles.size(); host_species_i++){
        std::unordered_map<int, int> allele_dist = hostPool.getAlleleDistributionAcrossAllLoci(host_species_i);
        hostAllelePool.purgeUnused(host_species_i, allele_dist);
    }

    for(int patho_species_i = 0; patho_species_i < pathogenAllelePool.alleles.size(); patho_species_i++){
        std::unordered_map<int, int> haplotype_dist = pathogenPool.getHaplotypeDistribution(patho_species_i);
        pathogenAllelePool.purgeUnused(patho_species_i, haplotype_dist);
    }

}

