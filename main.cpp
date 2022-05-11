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

    std::cout << "config: \n" << config.dump(4) << std::endl;

    auto pool_init_start = std::chrono::steady_clock::now();


    //TODO(JAN): smart pointer statt referenzen
    // config polymorphism -> instantiate object based on config
    // create shared ptr and pass that to the simulation env
    /*
    std::shared_ptr<InfectionRegime> infectionRegime;
    if(config["regime"] == "random"){
        infectionRegime = std::make_shared<RandomInfectionRegime>();
    }
    */

    RandomInfectionRegime randomInfectionRegime(config);
    SimulationEnvironment env(config, randomInfectionRegime);
    env.initialize();

    auto pool_init_end = std::chrono::steady_clock::now();

    std::cout << "pool init time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(pool_init_end - pool_init_start).count()
              << " ms" << std::endl;

    // sanity check
    int randomHostSpeciesIndex = rand() % env.hostPool.hosts.size();
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
    const int haplotype_id = randomPathogen.haplotypeId;
    const Allele& haplotype = env.pathogenAllelePool.alleles[randomPathogenSpeciesIndex][haplotype_id];

    std::cout << "haplotype verification " << haplotype.sequence << std::endl;

    const int allele_1_ht_merit = Helper::generate_merit(allele_1.sequence, haplotype.sequence);
    const int allele_1_ht_cached_merit = env.meritCache.get(randomHostSpeciesIndex, allele_1_id, randomPathogenSpeciesIndex, haplotype_id);

    std::cout << "allele 1 <-> haplotype merit: " << allele_1_ht_merit << std::endl;
    std::cout << "allele 1 <-> haplotype merit (cached): " << allele_1_ht_cached_merit << std::endl;

    auto simulation_start = std::chrono::steady_clock::now();

    env.step();

    auto simulation_end = std::chrono::steady_clock::now();

    randomHost.print();
    randomPathogen.print();


    std::cout << "Simulation Time: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(simulation_end - simulation_start).count()
         << " ms" << std::endl;

    return 0;
}

