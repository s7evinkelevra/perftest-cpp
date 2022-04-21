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
#include "InfectionRegiemes/InfectionRegieme.h"

#include "nlohmann/json.hpp"
using json = nlohmann::json;

class SimulationEnvironment {
private:
    InfectionRegieme& infectionRegieme;
public:
    SimulationEnvironment(json initialConfig, InfectionRegieme& infectionRegieme);
    json config;

    AllelePool hostAllelePool;
    AllelePool pathogenAllelePool;

    void initializeHostAllelePool();
    void initializePathogenAllelePool();

    MeritCache meritCache;

    void initializeMeritCache();

    HostPool hostPool;
    PathogenPool pathogenPool;

    void initializeHostPool();
    void initializePathogenPool();

    // stringify host
    void printHost(int index);
    void printPathogen(int index);


    // stat functions


};


#endif //PERFTEST_CPP_SIMULATIONENVIRONMENT_H
