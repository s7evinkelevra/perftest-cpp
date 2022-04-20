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
using json = nlohmann::json;

namespace po = boost::program_options;

#include "src/Host.h"
#include "src/HostPool.h"
#include "src/Pathogen.h"
#include "src/PathogenPool.h"
#include "src/Helper.h"
#include "src/AllelePool.h"
#include "src/Allele.h"
#include "src/SimulationEnvironment.h"


int main(int argc, char const *argv[]) {
    std::string configPath;

    // setup command line arguments
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("config-file,C", po::value<std::string>(&configPath), "path to the config file");

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

    SimulationEnvironment env(config);
    env.initializeHostAllelePool();
    env.initializePathogenAllelePool();

    env.initializeMeritCache();

    env.initializeHostPool();
    env.initializePathogenPool();

    // sanity check
    int randomHostIndex = rand() % env.hostPool.hosts.size();
    env.printHost(randomHostIndex);

    int randomPathogenIndex = rand() % env.pathogenPool.pathogens.size();
    env.printPathogen(randomPathogenIndex);


    auto pool_init_end = std::chrono::steady_clock::now();

    std::cout << "pool init time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(pool_init_end - pool_init_start).count()
              << " ms" << std::endl;

    auto simulation_start = std::chrono::steady_clock::now();

    for(int i = 0; i < 22500; i++) {
        //auto start_merit = std::chrono::steady_clock::now();
        const std::string allele_test = Helper::gen_random(9);
        const std::string haplotype_test = Helper::gen_random(2000);

        //auto result = edit_distance(allele_test, haplotype_test);

        const int merit = Helper::generate_merit(allele_test, haplotype_test);

        //auto end_merit = std::chrono::steady_clock::now();
        //std::cout << allele_test << " - " << haplotype_test << std::endl;
        //std::cout << "merit: " << merit << std::endl << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_merit - start_merit).count() << std::endl;
    }

    auto simulation_end = std::chrono::steady_clock::now();

    std::cout << "Simulation Time: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(simulation_end - simulation_start).count()
         << " ms" << std::endl;

    return 0;
}

