//
// Created by Jan on 18.04.2022.
//

#ifndef PERFTEST_CPP_PATHOGEN_H
#define PERFTEST_CPP_PATHOGEN_H

#include <string>

class Pathogen {
public:
    Pathogen(int pathogenId, double initialFitness, int speciesId, std::string pathogenHaplotype);
    int id;
    double fitness;
    int species;
    std::string haplotype;

    void print();
};


#endif //PERFTEST_CPP_PATHOGEN_H
