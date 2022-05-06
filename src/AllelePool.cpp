//
// Created by Jan on 18.04.2022.
//

#include "AllelePool.h"

#include <utility>
#include "Allele.h"

int AllelePool::addAllele(int species_id, std::string sequence) {
    unsigned long alleleId = alleles[species_id].size();
    alleles[species_id].push_back(Allele((int)alleleId,std::move(sequence)));
}
