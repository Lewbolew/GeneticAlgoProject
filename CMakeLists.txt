cmake_minimum_required(VERSION 3.7)
project(MyGA)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp Chromosome.cpp GeneticAlgo.cpp Log.cpp Population.cpp)
add_executable(MyGA ${SOURCE_FILES})
TARGET_LINK_LIBRARIES (MyGA pthread)
