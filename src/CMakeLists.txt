# GLOB_RECURSE will allow you to find everything with extension in subdirs

if (PROTOBUF_FOUND)
  file(GLOB sources *.cpp proto/*.cc)
  file(GLOB headers *.h proto/*.h)
else (PROTOBUF_FOUND)
	file(GLOB sources *.cpp)
	file(GLOB headers *.h)
endif (PROTOBUF_FOUND)

add_executable(express ${sources} ${headers})

target_link_libraries(express ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} "${CMAKE_CURRENT_SOURCE_DIR}/../bamtools/lib/libbamtools.a" "tcmalloc")

if(WIN32)
  target_link_libraries(express "${CMAKE_CURRENT_SOURCE_DIR}/../win_build/zlibd.lib")  
else(WIN32)  
  target_link_libraries(express "pthread")
endif(WIN32)

if (PROTOBUF_FOUND)
  target_link_libraries(express ${PROTOBUF_LIBRARIES})
endif(PROTOBUF_FOUND)

install(TARGETS express DESTINATION bin)
