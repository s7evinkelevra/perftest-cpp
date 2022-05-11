//
// Created by Jan on 11.05.2022.
//

#include "Random.h"

#include <random>
#include <chrono>

Random::Random() {
    // init random number generator
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
}

void Random::reseed(unsigned int seed) {
    rng.seed(seed);
}

unsigned int Random::sampleIntUniUnsignedInt(unsigned int from, unsigned int thru) {
    static std::uniform_int_distribution<unsigned int> d{};
    using parm_t = decltype(d)::param_type;
    return d( rng, parm_t{from, thru} );
}

float Random::sampleRealUniFloat(float from, float thru) {
    static std::uniform_real_distribution<float> d{};
    using parm_t = decltype(d)::param_type;
    return d( rng, parm_t{0, 1} );
}

double Random::sampleRealUniDouble(double from, double thru) {
    static std::uniform_real_distribution<double> dd{};
    using parm_t = decltype(dd)::param_type;
    return dd( rng, parm_t{from, thru} );
}

float Random::sampleGaussian(float mean, float variance) {
    static std::normal_distribution<float> d{};
    using parm_t = decltype(d)::param_type;
    return d( rng, parm_t{mean, variance} );
}
