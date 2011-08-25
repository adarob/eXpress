cmake_minimum_required (VERSION 2.6)

project(express)

set(${PROJECT_NAME}_VERSION_MAJOR 0)
set(${PROJECT_NAME}_VERSION_MINOR 1)
set(${PROJECT_NAME}_VERSION_PATCH 0)

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os ${CMAKE_CXX_FLAGS}")

set(CMAKE_BUILD_TYPE Release)

find_package(Boost 1.41
    COMPONENTS
    thread
	system
	filesystem
    program_options
    REQUIRED)

include_directories(${Boost_INCLUDE_DIRS} "bamtools/include")
link_directories ("../../bamtools/lib") 

add_subdirectory(src)
