# version requirements
cmake_minimum_required(VERSION 3.20.0)
cmake_policy(VERSION 3.20.0)

# which compiler does this actually use? ...
 message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")

# output
project(2MC VERSION 0.1.0 LANGUAGES CXX)
add_executable(2MC)
add_subdirectory(src)

# Link OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(2MC PRIVATE OpenMP::OpenMP_CXX)

# compilation options
set(CMAKE_CXX_STANDARD_REQUIRED "ON")
set_property(TARGET 2MC PROPERTY CXX_STANDARD 20)
