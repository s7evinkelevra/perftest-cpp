//
// Created by Jan on 18.04.2022.
//

#include "AllelePool.h"

#include <utility>
#include "Allele.h"

unsigned long AllelePool::addAllele(int species_id, int parent_id, std::string sequence) {
    unsigned long alleleId = alleles[species_id].size();
    alleles[species_id].push_back(Allele(parent_id, (int)alleleId,std::move(sequence)));
    return alleleId;
}

//TODO(JAN): completely untested
//BUG(JAN): ALLELES RELY ON INDEX==ID, THIS BREAKS THAT!!! DO NOT USE UNTIL SOLUTION IS FOUND
void AllelePool::purgeUnused(int species_id, std::unordered_map<int,int>& allele_dist){
    std::vector<Allele> used_alleles;
    used_alleles.reserve(allele_dist.size());

    for(auto& item : allele_dist){
        used_alleles.emplace_back(alleles[species_id][item.first]);
    }

    alleles[species_id] = used_alleles;
}