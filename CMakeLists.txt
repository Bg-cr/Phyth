cmake_minimum_required(VERSION 3.15)
project(Phyth)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(Phyth INTERFACE)

target_include_directories(Phyth INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

file(GLOB EXAMPLE_SOURCES "examples/*.cpp")

set(CMAKE_CONFIGURE_DEPENDS ${EXAMPLE_SOURCES})

foreach(EXAMPLE_SOURCE ${EXAMPLE_SOURCES})
    get_filename_component(TEST_NAME ${EXAMPLE_SOURCE} NAME_WE)
    add_executable(${TEST_NAME} ${EXAMPLE_SOURCE})
    target_link_libraries(${TEST_NAME} PRIVATE Phyth)
    target_include_directories(${TEST_NAME} PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
endforeach()