//
// Created by Jan on 20.04.2022.
//

#include "MeritCache.h"
#include <string>

std::string MeritCache::keyFromIds(int host_species_id, int host_allele_id, int patho_species_id, int patho_allele_id) {
    return std::to_string(host_species_id) + "_" + std::to_string(host_allele_id) + "_" + std::to_string(patho_species_id) + "_" + std::to_string(patho_allele_id);
}

int MeritCache::get(int host_species_id, int host_allele_id, int patho_species_id, int patho_allele_id) {
    std::string key = keyFromIds(host_species_id, host_allele_id, patho_species_id, patho_allele_id);
    return cache[key];
}

void MeritCache::set(int host_species_id, int host_allele_id, int patho_species_id, int patho_allele_id, int value) {
    std::string key = keyFromIds(host_species_id, host_allele_id, patho_species_id, patho_allele_id);
    cache[key] = value;
}


