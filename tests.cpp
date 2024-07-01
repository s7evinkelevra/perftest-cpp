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

void test_levenshtein(){
    std::string a = "teste";
    std::string b = "testa";

    std::cout << "levenshtein distance: " << Helper::LevenshteinDistance(a, b) << std::endl;
}

void test_hamming(){
    std::string a = "teste";
    std::string b = "testa";

    std::cout << "hamming distance: " << Helper::HammingDistance(a, b) << std::endl;
}

int main(int argc, char const *argv[]){
    unsigned int threads_available = std::thread::hardware_concurrency();
    std::cout << "threads available: " << threads_available << "\n";
    omp_set_num_threads((int)threads_available);

    test_levenshtein();
    test_hamming();

    const int N = 10000;

    long sum_lev = 0;

    for(int i = 0; i < N; i++) {
        auto haplotype = gen_random(2000);
        auto allele = gen_random(9);
        auto start_lev = std::chrono::steady_clock::now();

        Helper::generate_merit(allele, haplotype);
        auto end_lev = std::chrono::steady_clock::now();
        sum_lev += std::chrono::duration_cast<std::chrono::microseconds>(end_lev - start_lev).count();
    }


    long sum_hamming = 0;
    for(int i = 0; i < N; i++) {
        auto haplotype = gen_random(2000);
        auto allele = gen_random(9);

        auto start_hamming = std::chrono::steady_clock::now();
        Helper::generate_merit_hamming(allele, haplotype);
        auto end_hamming = std::chrono::steady_clock::now();
        sum_hamming += std::chrono::duration_cast<std::chrono::microseconds>(end_hamming - start_hamming).count();
    }



    std::cout << "lev run time: "
              << sum_lev
              << " us" << std::endl;

    std::cout << "hamming run time: "
              << sum_hamming
              << " us" << std::endl;

}