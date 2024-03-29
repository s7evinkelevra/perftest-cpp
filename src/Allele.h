//
// Created by Jan on 18.04.2022.
//

#ifndef PERFTEST_CPP_ALLELE_H
#define PERFTEST_CPP_ALLELE_H

#include <string>

class Allele {
public:
    Allele(int parent_id, int allele_id, int created_at_generation, std::string allele_sequence);
    int id;
    int parentId;
    std::string sequence;
    int createdAtGeneration;
};


#endif //PERFTEST_CPP_ALLELE_H
