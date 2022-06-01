#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <string_view>
#include <chrono>
#include <algorithm>
#include <vector>
#include <deque>
#include <boost/program_options.hpp>
#include <unordered_map>

#include "src/nlohmann/json.hpp"

#include "src/Host.h"
#include "src/HostPool.h"
#include "src/Pathogen.h"
#include "src/PathogenPool.h"
#include "src/Helper.h"
#include "src/AllelePool.h"
#include "src/Allele.h"
#include "src/SimulationEnvironment.h"
#include "src/InfectionRegiemes/RandomInfectionRegime.h"

using json = nlohmann::json;

namespace po = boost::program_options;




int main(int argc, char const *argv[]) {
    std::string configPath;

    // setup command line arguments
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("config-file,c", po::value<std::string>(&configPath), "path to the config file");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if(vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    if(vm.count("config-file")) {
        std::cout << "using config at path " << configPath << std::endl;
    }else{
        std::cout << "no config path set, using path config.json" << std::endl;
        configPath = "config.json";
    }

    // read config file
    std::ifstream configStream (configPath);
    if(!configStream){
        std::cout << "config file not found, exiting." << std::endl;
        return 1;
    }


    // parse config to json
    json config;
    configStream >> config;
    configStream.close();

    std::cout << "config: \n" << config.dump(4) << std::endl;

    auto pool_init_start = std::chrono::steady_clock::now();


    SimulationEnvironment env(config);
    env.initialize();

    auto pool_init_end = std::chrono::steady_clock::now();

    std::cout << "pool init time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(pool_init_end - pool_init_start).count()
              << " ms" << std::endl;

    // sanity check
    /* int randomHostSpeciesIndex = rand() % env.hostPool.hosts.size();
    int randomHostIndex = rand() % env.hostPool.hosts[randomHostSpeciesIndex].size();
    env.printHost(randomHostSpeciesIndex, randomHostIndex);

    Host& randomHost = env.hostPool.hosts[randomHostSpeciesIndex][randomHostIndex];
    int allele_1_id = randomHost.chromosome_1_allele_ids[0];
    int allele_2_id = randomHost.chromosome_2_allele_ids[0];

    const Allele& allele_1 = env.hostAllelePool.alleles[randomHostSpeciesIndex][allele_1_id];
    const Allele& allele_2 = env.hostAllelePool.alleles[randomHostSpeciesIndex][allele_2_id];

    std::cout << "allele 1 verification: " << allele_1.sequence << std::endl;
    std::cout << "allele 2 verification: " << allele_2.sequence << std::endl;

    int randomPathogenSpeciesIndex = rand() % env.pathogenPool.pathogens.size();
    int randomPathogenIndex = rand() % env.pathogenPool.pathogens[randomPathogenSpeciesIndex].size();
    env.printPathogen(randomPathogenSpeciesIndex, randomPathogenIndex);

    Pathogen& randomPathogen = env.pathogenPool.pathogens[randomPathogenSpeciesIndex][randomPathogenIndex];
    const int haplotype_id = randomPathogen.haplotype_id;
    const Allele& haplotype = env.pathogenAllelePool.alleles[randomPathogenSpeciesIndex][haplotype_id];

    std::cout << "haplotype verification " << haplotype.sequence << std::endl;

    const int allele_1_ht_merit = Helper::generate_merit(allele_1.sequence, haplotype.sequence);
    const int allele_1_ht_cached_merit = env.meritCache.get(randomHostSpeciesIndex, allele_1_id, randomPathogenSpeciesIndex, haplotype_id);

    std::cout << "allele 1 <-> haplotype merit: " << allele_1_ht_merit << std::endl;
    std::cout << "allele 1 <-> haplotype merit (cached): " << allele_1_ht_cached_merit << std::endl; */

    auto simulation_start = std::chrono::steady_clock::now();

    env.setBurnInMode();

    env.writeHostData();
    env.writeHostGenomeData();
    env.writeHostAlleleData();

    env.writePathogenData();
    env.writePathogenGenomeData();
    env.writePathogenAlleleData();

    env.writeMetaData();

    for(int burnin_generation = 0; burnin_generation < env.config["burnin_generations"]; burnin_generation++){
        env.step();
        if(burnin_generation % 10 == 0){
            env.purgeUnusedAlleles();
        }
    }

    env.writeHostData();
    env.writeHostGenomeData();
    env.writeHostAlleleData();

    env.writePathogenData();
    env.writePathogenGenomeData();
    env.writePathogenAlleleData();

    env.writeMetaData();

    std::cout << "allele distribution across all loci: \n";
    std::vector<std::unordered_map<int,int>> allele_dist_per_species = env.hostPool.getAlleleDistributionAcrossAllLociPerSpecies();

    for (int species_i = 0; species_i < allele_dist_per_species.size(); species_i++){
        std::cout << "---- species" << species_i << "----\n";
        for(auto& item : allele_dist_per_species[species_i]){
            std::cout << "allele " << item.first << " -> " << item.second << "\n";
        }
    }

    std::cout << "haplotype distribution \n";
    std::vector<std::unordered_map<int,int>> haplotype_dist_per_species = env.pathogenPool.getHaplotypeDistributionsPerSpecies();

    for (int species_i = 0; species_i < haplotype_dist_per_species.size(); species_i++){
        std::cout << "---- species" << species_i << "----\n";
        for(auto& item : haplotype_dist_per_species[species_i]){
            std::cout << "haplotype " << item.first << " -> " << item.second << "\n";
        }
    }

    std::cout << "before purge alleles in use in patho species 0: " << env.pathogenAllelePool.alleles[0].size() << "\n";
    env.purgeUnusedAlleles();
    std::cout << "after purge alleles in use in patho species 0: " << env.pathogenAllelePool.alleles[0].size() << "\n";

    env.setDefaultMode();
    for(int generation = 0; generation < env.config["generations"]; generation++){
        env.step();
        env.writeHostAlleleData();

        env.writePathogenAlleleData();

        env.writeMetaData();

        if(generation % 20 == 0){
            env.purgeUnusedAlleles();
            env.writeHostData();
            env.writeHostGenomeData();

            env.writePathogenData();
            env.writePathogenGenomeData();

        }
    }





    auto simulation_end = std::chrono::steady_clock::now();

    //randomHost.print();
    //randomPathogen.print();


    std::cout << "total simulation run time: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(simulation_end - simulation_start).count()
         << " ms" << std::endl;

    return 0;
}

