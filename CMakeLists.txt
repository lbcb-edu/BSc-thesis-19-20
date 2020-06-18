cmake_minimum_required(VERSION 3.10)

project(HERAalgorithm VERSION 1.0
            DESCRIPTION "BSc thesis 2020, HERA bioinformatics algorithm implementation"
            LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(sequences EXCLUDE_FROM_ALL)

add_executable(hera src/main.cpp src/preklapanje.cpp src/parsiranje.cpp src/gradi.cpp src/graf_a.cpp src/graf_im.cpp)

target_link_libraries(hera PUBLIC
                    bioparser
                    sequences)

target_include_directories(hera PUBLIC
                            sequences)