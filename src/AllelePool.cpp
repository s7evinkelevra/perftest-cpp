//
// Created by Jan on 18.04.2022.
//

#include "AllelePool.h"

#include <utility>
#include <iostream>
#include "Allele.h"

unsigned long AllelePool::addAllele(int species_id, int parent_id, int created_at_generation, std::string sequence) {
    int& allele_count = total_allele_counts_per_species[species_id];
    allele_count++;
    alleles[species_id].emplace(allele_count, Allele(parent_id, (int)allele_count, created_at_generation, std::move(sequence)));
    return allele_count;
}

void AllelePool::purgeUnused(int species_id, std::unordered_map<int,int>& allele_dist){
    std::unordered_map<int, Allele> used_alleles;
    used_alleles.reserve(allele_dist.size());

    //std::cout << "species id: " << species_id << "\n";
    //std::cout << "testing allele 10: " << alleles[species_id].at(10).id << "\n";

    for(auto& item : allele_dist){
        //std::cout << "allele: " << item.first << " count: " << item.second << "\n";
        Allele& used_allele = alleles[species_id].at(item.first);
        used_alleles.insert(std::make_pair(item.first, used_allele));
    }


    alleles[species_id] = used_alleles;
}