cmake_policy(SET CMP0048 NEW)
cmake_minimum_required(VERSION 3.0)

get_filename_component(filename generate_bprime_table NAME)

project(${filename})

# Let's require mutation++. It should export a suitable config that also
# exports include directories for Eigen. The `if` is needed to allow the
# examples to build during the process of building the entire mutation++
# tree.
if(NOT TARGET mutation++)
    find_package(mutation++ REQUIRED)
endif()

add_executable(${filename} ${filename}.cpp)
target_link_libraries(${filename} $ENV{MPP_DIRECTORY}/install/lib64/libmutation++.so)
