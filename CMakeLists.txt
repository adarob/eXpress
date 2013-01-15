cmake_minimum_required (VERSION 2.6)

project(express)

set(${PROJECT_NAME}_VERSION_MAJOR 1)
set(${PROJECT_NAME}_VERSION_MINOR 2)
set(${PROJECT_NAME}_VERSION_PATCH 3)

set(CMAKE_CXX_FLAGS "-Wall")

set(CMAKE_CXX_FLAGS_DEBUG "-g ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RHDEBINFO "-O3 -g ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os ${CMAKE_CXX_FLAGS}")

set(CMAKE_BUILD_TYPE Debug)
set(Boost_USE_STATIC_LIBS ON)

find_package(Boost 1.39
    COMPONENTS
    thread
	system
	filesystem
    program_options
    REQUIRED)

find_package(Protobuf)

if(WIN32)
	set(CMAKE_CXX_FLAGS "/EHsc")
	set(WIN32_INT 1)
else(WIN32)
	find_package(ZLIB REQUIRED)
	set(WIN32_INT 0)
endif(WIN32)

configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_SOURCE_DIR}/src/config.h"
  )

if (PROTOBUF_FOUND)
  include_directories(${Boost_INCLUDE_DIRS} ${PROTOBUF_INCLUDE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/bamtools/include") 
else (PROTOBUF_FOUND)
  include_directories(${Boost_INCLUDE_DIRS} "${CMAKE_CURRENT_SOURCE_DIR}/bamtools/include") 
endif(PROTOBUF_FOUND)
add_subdirectory(src)
