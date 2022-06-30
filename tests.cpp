#include <iostream>
#include <string>
#include <chrono>
#include <thread>
#include <omp.h>

#include <unordered_map>


#include "src/Helper.h"


std::string gen_random(const int len) {
    static const char AS[] =
            "ARNDCQGEHILKMFPSTWYV";
    std::string tmp_s;
    tmp_s.reserve(len);

    for (int i = 0; i < len; ++i) {
        tmp_s += AS[rand() % (sizeof(AS) - 1)];
    }

    return tmp_s;
}

int main(int argc, char const *argv[]){
    std::string a = "tester";
    std::string b = "testa";

    unsigned int threads_available = std::thread::hardware_concurrency();
    std::cout << "threads available: " << threads_available << "\n";
    omp_set_num_threads(threads_available);

    auto start_time = std::chrono::steady_clock::now();

    std::unordered_map<int,int> dist;

    int total = 100000;
    for(int i = 0; i < total; i++) {
        std::string random_peptide = gen_random(9);
        std::string random_haplotype = gen_random(2000);
        std::string random_mhc_mol = gen_random(9);

        //const int merit = Helper::LevenshteinDistance(random_mhc_mol, random_peptide);
        const int merit = Helper::generate_merit(random_mhc_mol, random_haplotype);

        dist[merit]++;
    }

    for(const auto& item : dist) {
        std::cout << "merit " << item.first << " -> " << item.second << " | " << (float)item.second / (float)total << std::endl;
    }

    auto end_time = std::chrono::steady_clock::now();

    std::cout << "run time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
              << " ms" << std::endl;

}