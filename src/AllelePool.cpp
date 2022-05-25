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
