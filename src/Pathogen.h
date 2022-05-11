//
// Created by Jan on 18.04.2022.
//

#ifndef PERFTEST_CPP_PATHOGEN_H
#define PERFTEST_CPP_PATHOGEN_H

#include <string>

class Pathogen {
public:
    Pathogen(int pathogenId, double initialFitness, int speciesId, int pathogenHaplotypeId);
    int id;
    double fitness;
    int species;
    int haplotypeId;

    int infection_count;
    int no_infection_count;

    void updateFitness();

    void print();
};


#endif //PERFTEST_CPP_PATHOGEN_H
