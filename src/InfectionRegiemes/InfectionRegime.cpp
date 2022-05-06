//
// Created by Jan on 21.04.2022.
//

#include "InfectionRegime.h"

#include <iostream>
#include <utility>


InfectionRegime::InfectionRegime(json initialConfig) {
    config = std::move(initialConfig);
}

void InfectionRegime::testMethod() {
    std::cout << "base class implementation" << std::endl;
}

void InfectionRegime::infect() {
    std::cout << "base class infect" << std::endl;
}
