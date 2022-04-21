cmake_minimum_required(VERSION 3.16)
project(perftest_cpp)

set(CMAKE_CXX_STANDARD 20)

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
set(Boost_USE_STATIC_LIBS   ON)

# copy config file to build directory for easy consumption during dev
configure_file(config.json config.json COPYONLY)
#set(CMAKE_CXX_FLAGS "-O3")
add_executable(perftest_cpp main.cpp src/HostPool.cpp src/HostPool.h src/Host.cpp src/Host.h src/Allele.cpp src/Allele.h src/Helper.h src/AllelePool.cpp src/AllelePool.h src/Pathogen.cpp src/Pathogen.h src/PathogenPool.cpp src/PathogenPool.h src/SimulationEnvironment.cpp src/SimulationEnvironment.h src/MeritCache.cpp src/MeritCache.h src/Helper.cpp src/InfectionRegiemes/InfectionRegieme.cpp src/InfectionRegiemes/InfectionRegieme.h src/InfectionRegiemes/RandomInfectionRegieme.cpp src/InfectionRegiemes/RandomInfectionRegieme.h)

TARGET_LINK_LIBRARIES(perftest_cpp ${Boost_LIBRARIES})