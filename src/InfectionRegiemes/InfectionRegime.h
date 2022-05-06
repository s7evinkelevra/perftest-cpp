//
// Created by Jan on 21.04.2022.
//

#ifndef PERFTEST_CPP_INFECTIONREGIME_H
#define PERFTEST_CPP_INFECTIONREGIME_H

#include "../nlohmann/json.hpp"
using json = nlohmann::json;

class InfectionRegime {
public:
    InfectionRegime(json initialConfig);

    json config;

    virtual void testMethod();
    virtual void infect() = 0;
};


#endif //PERFTEST_CPP_INFECTIONREGIME_H
