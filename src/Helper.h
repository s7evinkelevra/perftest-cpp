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

    template<typename T> inline
    typename T::size_type HammingDistance(const T &source, const T &target) {
        // won't happen
        //if (source.size() != target.size()) {
        //    throw std::invalid_argument("Hamming distance is only defined for strings of equal length");
        //}

        using TSizeType = typename T::size_type;
        TSizeType distance = 0;
        for (TSizeType i = 0; i < source.size(); ++i) {
            if (source[i] != target[i]) {
                ++distance;
            }
        }

        return distance;
    }

    std::string gen_random(const int len);

    int generate_merit(std::string_view allele, std::string_view haplotype);
    int generate_merit_hamming(std::string_view allele, std::string_view hyplotype);

    std::vector<std::string> split(std::string_view s, char delimiter);

}


#endif //PERFTEST_CPP_HELPER_H
