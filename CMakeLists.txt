cmake_minimum_required(VERSION 3.5)

set (CMAKE_CXX_STANDARD 11)

project(MDI_metadynamics
        LANGUAGES C CXX)

add_subdirectory(lib)
add_subdirectory(MDI_metadynamics)
