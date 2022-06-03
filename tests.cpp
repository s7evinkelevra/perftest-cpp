#include <iostream>
#include <string>
#include "src/Helper.h"

int main(int argc, char const *argv[]){
  std::string a = "tester";
  std::string b = "testa";

  for(int i = 0; i < 100000000; i++) {
      Helper::LevenshteinDistance(a,b);
  }


}