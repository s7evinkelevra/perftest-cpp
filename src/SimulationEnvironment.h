//
// Created by Jan on 20.04.2022.
//

#ifndef PERFTEST_CPP_SIMULATIONENVIRONMENT_H
#define PERFTEST_CPP_SIMULATIONENVIRONMENT_H

#include <random>

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
    InfectionRegime& infectionRegieme;

    void initializeHostAllelePool();
    void initializePathogenAllelePool();
    void initializeMeritCache();
    void initializeHostPool();
    void initializePathogenPool();

public:
    SimulationEnvironment(json initialConfig, InfectionRegime& infectionRegieme);
    json config;
    //TODO(JAN): dont forget to use these member vars
    int generation;
    int pathogenGeneration;
    bool bInfection;
    bool bHostMutation;
    bool bPathogenMutation;
    bool bReproduction;

    AllelePool hostAllelePool;
    AllelePool pathogenAllelePool;

    //TODO(JAN): control access through getter/setter
    MeritCache meritCache;

    HostPool hostPool;
    PathogenPool pathogenPool;


    void initialize();

    // stringify host
    void printHost(int species, int index);
    void printPathogen(int species, int index);

    void testMethod();

    // stat functions


    // simulation
    void step();

};


#endif //PERFTEST_CPP_SIMULATIONENVIRONMENT_H
