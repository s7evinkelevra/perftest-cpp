//
// Created by Jan on 18.04.2022.
//

#include "Allele.h"

#include <utility>
#include <chrono>


Allele::Allele(int parent_id, int alleleId, int created_at_generation, std::string allele_sequence) {
    id = alleleId;
    parentId = parent_id;
    sequence = std::move(allele_sequence);
    createdAtGeneration = created_at_generation;
}
