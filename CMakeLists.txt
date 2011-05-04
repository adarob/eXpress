cmake_minimum_required (VERSION 2.6)

project(cuffexpress)

set(${PROJECT_NAME}_VERSION_MAJOR 0)
set(${PROJECT_NAME}_VERSION_MINOR 01)

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os ${CMAKE_CXX_FLAGS}")

set(CMAKE_BUILD_TYPE Release)

find_package(Boost 1.42
    COMPONENTS
        thread
    REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(src)
