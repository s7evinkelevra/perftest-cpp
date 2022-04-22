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
public:
    SimulationEnvironment(json initialConfig, InfectionRegime& infectionRegieme);
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

    void testMethod();

    // stat functions


};


#endif //PERFTEST_CPP_SIMULATIONENVIRONMENT_H
