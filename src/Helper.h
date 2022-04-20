//
// Created by Jan on 18.04.2022.
//

#ifndef PERFTEST_CPP_HELPER_H
#define PERFTEST_CPP_HELPER_H


#include <string>
#include <vector>


namespace Helper{
    template<typename T> inline
    typename T::size_type LevenshteinDistance(const T &source, const T &target) {
        if (source.size() > target.size()) {
            return LevenshteinDistance(target, source);
        }

        using TSizeType = typename T::size_type;
        const TSizeType min_size = source.size(), max_size = target.size();
        std::vector<TSizeType> lev_dist(min_size + 1);

        for (TSizeType i = 0; i <= min_size; ++i) {
            lev_dist[i] = i;
        }

        for (TSizeType j = 1; j <= max_size; ++j) {
            TSizeType previous_diagonal = lev_dist[0], previous_diagonal_save;
            ++lev_dist[0];

            for (TSizeType i = 1; i <= min_size; ++i) {
                previous_diagonal_save = lev_dist[i];
                if (source[i - 1] == target[j - 1]) {
                    lev_dist[i] = previous_diagonal;
                } else {
                    lev_dist[i] = std::min(std::min(lev_dist[i - 1], lev_dist[i]), previous_diagonal) + 1;
                }
                previous_diagonal = previous_diagonal_save;
            }
        }

        return lev_dist[min_size];
    }

    std::string gen_random(const int len);

    int generate_merit(std::string_view allele, std::string_view haplotype);

}


#endif //PERFTEST_CPP_HELPER_H
