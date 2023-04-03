cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)
project(ntHash_testing)

############################################################
# Create a library
############################################################

#Generate the static library from the library sources
add_library(ntHash STATIC 
    src/ntHash/nthash.cpp src/ntHash/nthash_consts.cpp src/ntHash/nthash_lowlevel.cpp src/ntHash/status.cpp
)

target_include_directories(ntHash
    PUBLIC 
        ${PROJECT_SOURCE_DIR}/include/ntHash
)


############################################################
# Create an executable
############################################################

# Add an executable with the above sources

add_executable(test 
    ntHash_testing.cpp
)

# link the new hello_library target with the hello_binary target

target_link_libraries( test
    PRIVATE 
        ntHash
)