cmake_minimum_required(VERSION 3.15)
project(perftest_cpp)

set(CMAKE_CXX_STANDARD 17)


find_package(OpenMP REQUIRED)

# copy config file to build directory for easy consumption during dev
configure_file(config.json config.json COPYONLY)
#set(CMAKE_CXX_FLAGS "-O3")
add_executable(perftest_cpp main.cpp src/HostPool.cpp src/HostPool.h src/Host.cpp src/Host.h src/Allele.cpp src/Allele.h src/Helper.h src/AllelePool.cpp src/AllelePool.h src/Pathogen.cpp src/Pathogen.h src/PathogenPool.cpp src/PathogenPool.h src/SimulationEnvironment.cpp src/SimulationEnvironment.h src/MeritCache.cpp src/MeritCache.h src/Helper.cpp src/InfectionRegiemes/InfectionRegime.cpp src/InfectionRegiemes/InfectionRegime.h src/InfectionRegiemes/RandomInfectionRegime.cpp src/InfectionRegiemes/RandomInfectionRegime.h src/Random.cpp src/Random.h src/CSVWriter.cpp src/CSVWriter.h)
add_executable(testing_cpp tests.cpp src/HostPool.cpp src/HostPool.h src/Host.cpp src/Host.h src/Allele.cpp src/Allele.h src/Helper.h src/AllelePool.cpp src/AllelePool.h src/Pathogen.cpp src/Pathogen.h src/PathogenPool.cpp src/PathogenPool.h src/SimulationEnvironment.cpp src/SimulationEnvironment.h src/MeritCache.cpp src/MeritCache.h src/Helper.cpp src/InfectionRegiemes/InfectionRegime.cpp src/InfectionRegiemes/InfectionRegime.h src/InfectionRegiemes/RandomInfectionRegime.cpp src/InfectionRegiemes/RandomInfectionRegime.h src/Random.cpp src/Random.h src/CSVWriter.cpp src/CSVWriter.h)

TARGET_LINK_LIBRARIES(perftest_cpp OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(testing_cpp OpenMP::OpenMP_CXX)
