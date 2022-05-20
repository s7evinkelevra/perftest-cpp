//
// Created by Jan on 11.05.2022.
//

#ifndef PERFTEST_CPP_RANDOM_H
#define PERFTEST_CPP_RANDOM_H

#include <random>



class Random {
public:
    std::mt19937_64 rng;

    Random();
    void reseed(unsigned int seed);

    // sample uniform dist
    unsigned int sampleIntUniUnsignedInt(unsigned int from, unsigned int thru);
    float sampleRealUniFloat(float from, float thru);
    double sampleRealUniDouble(double from, double thru);

    // sample gaussian
    float sampleGaussian(float mean, float variance);

    // sample binomial
    int sampleBinomial(int upper, double probability);



};


#endif //PERFTEST_CPP_RANDOM_H
