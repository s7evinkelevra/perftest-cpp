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
#include "CSVWriter.h"
#include <memory>


#include "nlohmann/json.hpp"
using json = nlohmann::json;

class SimulationEnvironment {
private:
    Random rng;


    std::unique_ptr<CSVWriter> hostDataCSV;
    std::unique_ptr<CSVWriter> hostGenomeDataCSV;
    std::unique_ptr<CSVWriter> hostAlleleDataCSV;
    std::unique_ptr<CSVWriter> hostLocusDataCSV;

    std::unique_ptr<CSVWriter> pathogenDataCSV;
    std::unique_ptr<CSVWriter> pathogenGenomeDataCSV;
    std::unique_ptr<CSVWriter> pathogenAlleleDataCSV;
    std::unique_ptr<CSVWriter> pathogenLocusDataCSV;


    void initializeHostAllelePool();
    void initializePathogenAllelePool();
    void initializeMeritCache();
    void initializeHostPool();
    void initializePathogenPool();

    void initializeCSVFiles();

public:
    SimulationEnvironment(json initialConfig);
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

    // stringify host
    void printHost(int species, int index);
    void printPathogen(int species, int index);

    // stat  functions

    // writing data
    void writeHostData();
    void writeHostGenomeData();
    void writeHostAlleleData();
    void writeHostLocusData();

    void writePathogenData();
    void writePathogenGenomeData();
    void writePathogenAlleleData();
    void writePathogenLocusData();

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
