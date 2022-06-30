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
#include <chrono>


#include "nlohmann/json.hpp"
using json = nlohmann::json;

class SimulationEnvironment {
private:
    Random rng;

    std::chrono::time_point<std::chrono::steady_clock> lastStepStart;
    std::chrono::time_point<std::chrono::steady_clock> lastStepEnd;


    std::unique_ptr<CSVWriter> hostDataCSV;
    std::unique_ptr<CSVWriter> hostGenomeDataCSV;
    std::unique_ptr<CSVWriter> hostAlleleDataCSV;
    std::unique_ptr<CSVWriter> hostAlleleSequenceDataCSV;
    std::unique_ptr<CSVWriter> hostLocusDataCSV;

    std::unique_ptr<CSVWriter> pathogenDataCSV;
    std::unique_ptr<CSVWriter> pathogenGenomeDataCSV;
    std::unique_ptr<CSVWriter> pathogenAlleleDataCSV;
    std::unique_ptr<CSVWriter> pathogenAlleleSequenceDataCSV;
    std::unique_ptr<CSVWriter> pathogenLocusDataCSV;

    std::unique_ptr<CSVWriter> metaDataCSV;


    void initializeHostAllelePool();
    void initializePathogenAllelePool();
    void initializeMeritCache();
    void initializeHostPool();
    void initializePathogenPool();

    void initializeOutputFiles();

public:
    SimulationEnvironment(json initialConfig);
    json config;

    unsigned int thread_count;

    int totalHostGenerations;
    int totalPathogenGenerations;

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
    void writeHostAlleleSequenceData();
    void writeHostLocusData();

    void writePathogenData();
    void writePathogenGenomeData();
    void writePathogenAlleleData();
    void writePathogenAlleleSequenceData();
    void writePathogenLocusData();

    void writeMetaData();

    void writeAllData();

    void flushAllDataToDisk();
    // housekeeping
    void purgeUnusedAlleles();

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
