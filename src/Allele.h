//
// Created by Jan on 18.04.2022.
//

#ifndef PERFTEST_CPP_ALLELE_H
#define PERFTEST_CPP_ALLELE_H

#include <string>

class Allele {
public:
    Allele(int alleleId, std::string allele_sequence);
    int id;
    std::string sequence;
    long createdAt;
};


#endif //PERFTEST_CPP_ALLELE_H
