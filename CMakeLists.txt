# version requirements
cmake_minimum_required(VERSION 3.20.0)
cmake_policy(VERSION 3.20.0)

# output
project(2MC VERSION 0.1.0 LANGUAGES CXX)
add_executable(2MC main.cpp)

# compilation options
set(CMAKE_CXX_STANDARD_REQUIRED "ON")
set_property(TARGET 2MC PROPERTY CXX_STANDARD 20)

# which compiler does this actually use? ...