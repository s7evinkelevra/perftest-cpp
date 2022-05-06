//
// Created by Jan on 21.04.2022.
//

#ifndef PERFTEST_CPP_RANDOMINFECTIONREGIME_H
#define PERFTEST_CPP_RANDOMINFECTIONREGIME_H


#include "InfectionRegime.h"

#include <utility>
#include <iostream>

class RandomInfectionRegime : public InfectionRegime {
public:
    RandomInfectionRegime(json initialConfig) : InfectionRegime(std::move(initialConfig)) {
        testInt = 0;
    };
    int testInt;
    void testMethod() override;
    void infect() override;
};


#endif //PERFTEST_CPP_RANDOMINFECTIONREGIME_H
