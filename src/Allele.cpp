//
// Created by Jan on 18.04.2022.
//

#include "Allele.h"

#include <utility>
#include <chrono>


Allele::Allele(int parent_id, int alleleId, std::string allele_sequence) {
    id = alleleId;
    parentId = parent_id;
    sequence = std::move(allele_sequence);
    const auto now = std::chrono::system_clock::now();
    createdAt = std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()).count();
}
