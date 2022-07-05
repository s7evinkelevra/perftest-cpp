//
// Created by Jan on 18.04.2022.
//

#ifndef PERFTEST_CPP_PATHOGEN_H
#define PERFTEST_CPP_PATHOGEN_H

#include <string>

class Pathogen {
public:
    Pathogen(int parentId, int pathogenId, double fitnessMinimum, int speciesId, int pathogenHaplotypeId);
    int parent_id;
    int id;
    double fitness;
    double fitness_minimum;
    int species;
    int haplotype_id;

    int infection_count;
    int no_infection_count;

    void updateFitness();
    void resetInfectionCount(int initialFitness);

    void print();
};


#endif //PERFTEST_CPP_PATHOGEN_H
