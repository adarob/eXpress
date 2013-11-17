cmake_minimum_required (VERSION 2.6)

project(express)

set(${PROJECT_NAME}_VERSION_MAJOR 1)
set(${PROJECT_NAME}_VERSION_MINOR 4)
set(${PROJECT_NAME}_VERSION_PATCH 4)

set(CMAKE_CXX_FLAGS "-Wall")

set(CMAKE_CXX_FLAGS_DEBUG "-g ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RHDEBINFO "-O3 -g ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os ${CMAKE_CXX_FLAGS}")

set(CMAKE_BUILD_TYPE Release)
set(Boost_USE_STATIC_LIBS ON)

find_package(Boost 1.39
    COMPONENTS
        thread
   	system
    	filesystem
    	program_options
    	date_time
    REQUIRED)

find_library(GPERFTOOLS_TCMALLOC_LIB tcmalloc)
if (GPERFTOOLS_TCMALLOC_LIB)
   message(STATUS "Found GPERFTOOLS_TCMALLOC: ${GPERFTOOLS_TCMALLOC_LIB}")
else (GPERFTOOLS_TCMALLOC_LIB)
   message(STATUS "Could NOT find GPERFTOOLS_TCMALLOC: Install to improve speed.")
endif(GPERFTOOLS_TCMALLOC_LIB)

find_package(Protobuf)
if (PROTOBUF_FOUND)
  include_directories(${Boost_INCLUDE_DIRS} ${PROTOBUF_INCLUDE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/bamtools/include") 
	set(PROTO_INT 1)
else (PROTOBUF_FOUND)
  include_directories(${Boost_INCLUDE_DIRS} "${CMAKE_CURRENT_SOURCE_DIR}/bamtools/include") 
	set(PROTO_INT 0)
endif(PROTOBUF_FOUND)

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


add_subdirectory(src)
