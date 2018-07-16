file(GLOB sources *.cpp)

if (PROTOBUF_FOUND) 
  file(GLOB ProtoFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.proto")
  PROTOBUF_GENERATE_CPP(PROTO_SOURCES PROTO_HEADERS ${ProtoFiles})
endif(PROTOBUF_FOUND)

add_executable(express ${sources} ${headers} ${PROTO_SOURCES} ${PROTO_HEADERS})

set(LIBRARIES ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${BAMTOOLS_LIBRARIES})

if (GPERFTOOLS_TCMALLOC) 
   set(LIBRARIES ${LIBRARIES} "libtcmalloc_minimal.a")
endif (GPERFTOOLS_TCMALLOC)

if (PROTOBUF_FOUND)
  get_filename_component(PROTOBUF_LIB_DIR ${PROTOBUF_LIBRARY} DIRECTORY)
  set(LIBRARIES ${LIBRARIES} "${PROTOBUF_LIB_DIR}/libprotobuf.a")
endif(PROTOBUF_FOUND)

if (WIN32)
target_link_libraries(express ${LIBRARIES})
elseif(APPLE)
target_link_libraries(express ${LIBRARIES} pthread)
elseif(UNIX)
target_link_libraries(express ${LIBRARIES} pthread rt)
endif()
install(TARGETS express DESTINATION bin)
