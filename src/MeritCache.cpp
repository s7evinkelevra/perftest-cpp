//
// Created by Jan on 20.04.2022.
//

#include "MeritCache.h"
#include <string>
#include <iostream>

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

float MeritCache::getAverage() {
    unsigned long size = cache.size();
    int sum = 0;

    for(auto& item : cache){
        sum += item.second;
    }

    return (float)sum/(float)size;
}

std::unordered_map<int, int> MeritCache::getDistribution() {
    std::unordered_map<int, int> dist;

    for(auto& item : cache){
        dist[item.second]++;
    }

/*    for(auto& item: dist){
        std::cout << item.first << " : " << item.second << std::endl;
    }*/

    return dist;
}

int MeritCache::getLowest() {
    int lowest = 9999;

    for(auto& item : cache){
        if(item.second < lowest){
            lowest = item.second;
        }
    }

    return lowest;
}

int MeritCache::getHighest() {
    int highest = -1;

    for(auto& item : cache){
        if(item.second > highest){
            highest = item.second;
        }
    }
    return highest;
}

bool MeritCache::exists(int host_species_id, int host_allele_id, int patho_species_id, int patho_allele_id) {
    std::string key = keyFromIds(host_species_id, host_allele_id, patho_species_id, patho_allele_id);
    return cache.find(key) != cache.end();
}


