//
// Created by Jan on 20.04.2022.
//

#ifndef PERFTEST_CPP_SIMULATIONENVIRONMENT_H
#define PERFTEST_CPP_SIMULATIONENVIRONMENT_H

#include "Random.h"

#include "HostPool.h"
#include "PathogenPool.h"
#include "AllelePool.h"
#include "Helper.h"
#include "MeritCache.h"
#include "InfectionRegiemes/InfectionRegime.h"

#include "nlohmann/json.hpp"
using json = nlohmann::json;

class SimulationEnvironment {
private:
    Random rng;
    InfectionRegime& infectionRegieme;

    void initializeHostAllelePool();
    void initializePathogenAllelePool();
    void initializeMeritCache();
    void initializeHostPool();
    void initializePathogenPool();

public:
    SimulationEnvironment(json initialConfig, InfectionRegime& infectionRegieme);
    json config;
    int totalHostGenerations;
    int totalPathogenGenerations;

    //TODO(JAN): dont forget to use these member vars
    bool bInfection;
    bool bHostFitnessproportionalReproduction;
    bool bPathogenFitnessproportionalReproduction;
    bool bHostMutation;
    bool bPathogenMutation;

    AllelePool hostAllelePool;
    AllelePool pathogenAllelePool;

    //TODO(JAN): control access through getter/setter
    MeritCache meritCache;

    HostPool hostPool;
    PathogenPool pathogenPool;


    void initialize();
    std::string generateSequence(int length);
    int addHostAllele(int host_species_id, const std::string& sequence);
    int addPathogenAllele(int patho_species_id, const std::string& sequence);

    // stringify host
    void printHost(int species, int index);
    void printPathogen(int species, int index);

    void testMethod();

    // stat functions


    // simulation
    void setDefaultMode();
    void setBurnInMode();
    void setNoCoevolutionMode();

    void step();

    void pathogenGeneration();

    void hostReproduction();

    void hostReproductionRandom();

    void hostMutation();

    void hostGeneration();

    void infection();

    void pathogenReproduction();

    void pathogenReproductionRandom();

    void pathogenMutation();
};


#endif //PERFTEST_CPP_SIMULATIONENVIRONMENT_H
