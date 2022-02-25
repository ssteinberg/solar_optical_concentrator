cmake_minimum_required(VERSION 3.20.0)
cmake_policy(VERSION 3.20.0)

# Output
project(tracer VERSION 0.1.0 LANGUAGES CXX)
add_executable(tracer src/main.cpp src/raytracer.cpp src/random_generator.cpp)

# Compilation options
set(CMAKE_CXX_STANDARD_REQUIRED "ON")
set_property(TARGET tracer PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
set_property(TARGET tracer PROPERTY CXX_STANDARD 20)
option (FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." True)

if (${FORCE_COLORED_OUTPUT})
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
       add_compile_options (-fdiagnostics-color=always)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
       add_compile_options (-fcolor-diagnostics)
    endif ()
endif ()

if (WIN32)
	add_definitions(-D_WIN32)
endif (WIN32)
if (UNIX)
	target_compile_options(tracer PRIVATE 
		-march=native -pipe 
		-mavx2
	)
    target_compile_options(tracer PUBLIC "$<$<CONFIG:RELEASE>:-O3>")

	add_definitions(-D_linux -D__AVX2__)

	target_link_libraries(tracer PUBLIC
		stdc++fs
		pthread atomic
	)
endif (UNIX)

# Include paths
target_include_directories(tracer PRIVATE 
    include/
    external/glm
)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

