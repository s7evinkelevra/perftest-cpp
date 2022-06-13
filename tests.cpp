#include <iostream>
#include <string>
#include <chrono>
#include <thread>
#include <omp.h>
#include "src/Helper.h"

int main(int argc, char const *argv[]){
    std::string a = "tester";
    std::string b = "testa";

    unsigned int threads_available = std::thread::hardware_concurrency();
    std::cout << "threads available: " << threads_available << "\n";
    omp_set_num_threads(threads_available);


    auto start_time = std::chrono::steady_clock::now();
    #pragma omp parallel for default(none) shared(a, b)
    for(int i = 0; i < 100000000; i++) {
        Helper::LevenshteinDistance(a,b);
    }

    auto end_time = std::chrono::steady_clock::now();

    std::cout << "run time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
              << " ms" << std::endl;

}