#include <iostream>
#include <string>
#include <chrono>
#include "src/Helper.h"

int main(int argc, char const *argv[]){
    std::string a = "tester";
    std::string b = "testa";

    auto start_time = std::chrono::steady_clock::now();

    for(int i = 0; i < 100000000; i++) {
        Helper::LevenshteinDistance(a,b);
    }

    auto end_time = std::chrono::steady_clock::now();

    std::cout << "run time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
              << " ms" << std::endl;

}