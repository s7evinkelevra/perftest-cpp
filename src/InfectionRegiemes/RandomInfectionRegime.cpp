//
// Created by Jan on 21.04.2022.
//

#include "RandomInfectionRegime.h"

#include <iostream>

void RandomInfectionRegime::testMethod() {
    std::cout << "derived class implementation, testint: " << testInt << std::endl;
}

void RandomInfectionRegime::infect() {
    std::cout << "testing config access (host count): " << config["hosts"]["n"] << std::endl;
    std::cout << "random infection regime infecting" << std::endl;
}
