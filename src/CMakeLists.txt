# GLOB_RECURSE will allow you to find everything with extension in subdirs
file(GLOB sources *.cpp)
file(GLOB headers *.h)

file(GLOB to_remove express-spark.cpp)
list(REMOVE_ITEM sources ${to_remove})

add_executable(express ${sources} ${headers})

if(WIN32)
  target_link_libraries(express ${Boost_LIBRARIES} "${CMAKE_CURRENT_SOURCE_DIR}/../win_build/zlibd.lib" "${CMAKE_CURRENT_SOURCE_DIR}/../bamtools/lib/libbamtools.lib")  
else(WIN32)  
  target_link_libraries(express ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} "${CMAKE_CURRENT_SOURCE_DIR}/../bamtools/lib/libbamtools.a" "pthread")
endif(WIN32)

install(TARGETS express DESTINATION bin)
