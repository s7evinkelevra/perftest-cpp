//
// Created by Jan on 20.04.2022.
//

#ifndef PERFTEST_CPP_MERITCACHE_H
#define PERFTEST_CPP_MERITCACHE_H

#include <deque>
#include <vector>

#include <unordered_map>
#include <string>

class MeritCache {
public:
    std::unordered_map<std::string, int> cache;
    void set(int host_species_id, int host_allele_id, int patho_species_id, int patho_allele_id, int value);
    int get(int host_species_id, int host_allele_id, int patho_species_id, int patho_allele_id);
    bool exists(int host_species_id, int host_allele_id, int patho_species_id, int patho_allele_id);

    float getAverage();
    std::unordered_map<int, int> getDistribution();
    int getLowest();
    int getHighest();

    static std::string keyFromIds(int host_species_id, int host_allele_id, int patho_species_id, int patho_haplotype_id);
};


#endif //PERFTEST_CPP_MERITCACHE_H
