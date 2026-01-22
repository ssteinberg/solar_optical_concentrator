add_subdirectory(geometry)

target_sources(2MC
        PRIVATE
        main.cpp
        bounding_box.cpp
        design.cpp
        node.cpp
        path.cpp
        tree.cpp
)

target_include_directories(2MC
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include)
