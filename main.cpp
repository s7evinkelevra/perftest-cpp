#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <chrono>

#include <algorithm>
#include <vector>

template<typename T>
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

/*
static int edit_distance(std::string_view s1, std::string_view s2) {

    const std::size_t len1 = s1.size();
    const std::size_t len2 = s2.size();

    std::vector<std::vector<int>> d(len1 + 1, std::vector<int>(len2 + 1));


    d[0][0] = 0;
    for(int i = 1; i <= len1; ++i) d[i][0] = i;
    for(int i = 1; i <= len2; ++i) d[0][i] = i;

    for(int i = 1; i <= len1; ++i)
        for(int j = 1; j <= len2; ++j)
            // note that std::min({arg1, arg2, arg3}) works only in C++11,
            // for C++98 use std::min(std::min(arg1, arg2), arg3)
            d[i][j] = std::min({d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1)});
    return d[len1][len2];
}
*/

std::string gen_random(const int len) {
    static const char AS[] =
            "ACDEFGHIKLMNOPQRSTUVWY";
    std::string tmp_s;
    tmp_s.reserve(len);

    for (int i = 0; i < len; ++i) {
        tmp_s += AS[rand() % (sizeof(AS) - 1)];
    }

    return tmp_s;
}

int generate_merit(std::string_view allele, std::string_view haplotype){
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

    return 1 - lowest_edit_distance / 9;
}



int main() {
    auto start = std::chrono::steady_clock::now();

    for(int i = 0; i < 22500; i++) {
        //auto start_merit = std::chrono::steady_clock::now();
        const std::string allele_test = gen_random(9);
        const std::string haplotype_test = gen_random(2000);

        //auto result = edit_distance(allele_test, haplotype_test);

        const int merit = generate_merit(allele_test, haplotype_test);

        //auto end_merit = std::chrono::steady_clock::now();
        //std::cout << allele_test << " - " << haplotype_test << std::endl;
        //std::cout << "merit: " << merit << std::endl << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_merit - start_merit).count() << std::endl;
    }

    auto end = std::chrono::steady_clock::now();

    std::cout << "Elapsed time in milliseconds: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
         << " ms" << std::endl;

    return 0;
}

