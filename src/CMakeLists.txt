if (PROTOBUF_FOUND)
  file(GLOB sources *.cpp proto/*.cc)
  file(GLOB headers *.h proto/*.h)
else (PROTOBUF_FOUND)
	file(GLOB sources *.cpp)
endif (PROTOBUF_FOUND)

add_executable(express ${sources} ${headers})

set(LIBRARIES ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})

if (GPERFTOOLS_TCMALLOC) 
   set(LIBRARIES ${LIBRARIES} "libtcmalloc_minimal.a")
endif (GPERFTOOLS_TCMALLOC)

if(WIN32)
  set(LIBRARIES ${LIBRARIES} "${CMAKE_CURRENT_SOURCE_DIR}/../bamtools/lib/libbamtools.lib" "${CMAKE_CURRENT_SOURCE_DIR}/../win_build/zlibd.lib")
else(WIN32)
  set(LIBRARIES ${LIBRARIES} "${CMAKE_CURRENT_SOURCE_DIR}/../bamtools/lib/libbamtools.a" "pthread")
endif(WIN32)

if (PROTOBUF_FOUND)
  set(LIBRARIES ${LIBRARIES} "libprotobuf.a")
endif(PROTOBUF_FOUND)

target_link_libraries(express ${LIBRARIES})

install(TARGETS express DESTINATION bin)
