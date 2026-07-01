cmake_minimum_required(VERSION 3.16)
project(QEC C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# QuEST as proper dependency
add_subdirectory(external/QuEST/QuEST)

# library layer 
add_library(qec_core
    src/encoder_functions.c
)

target_include_directories(qec_core PUBLIC include)

target_link_libraries(qec_core
    PUBLIC QuEST
    PUBLIC lapack
    PUBLIC lapacke
    PUBLIC blas
    PUBLIC m
)

# Final executable (experiment runner)
add_executable(qec_sim src/main.c)

target_link_libraries(qec_sim PRIVATE qec_core)
