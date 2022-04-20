//
// Created by Jan on 20.04.2022.
//

#include "Helper.h"


std::string Helper::gen_random(const int len) {
    static const char AS[] =
            "ACDEFGHIKLMNOPQRSTUVWY";
    std::string tmp_s;
    tmp_s.reserve(len);

    for (int i = 0; i < len; ++i) {
        tmp_s += AS[rand() % (sizeof(AS) - 1)];
    }

    return tmp_s;
}

int Helper::generate_merit(std::string_view allele, std::string_view haplotype){
    const long step_size = 1;
    const unsigned long window_width = allele.length();

    int lowest_edit_distance = -1;

    for(int i = 0; i < haplotype.length() - window_width + step_size; i += step_size) {
        std::string_view window = haplotype.substr(i,window_width);
        const int current_edit_distance = LevenshteinDistance(allele, window);
        if(lowest_edit_distance == -1 || lowest_edit_distance > current_edit_distance) {
            lowest_edit_distance = current_edit_distance;
        }
    }

    return lowest_edit_distance;
}