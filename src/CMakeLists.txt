# GLOB_RECURSE will allow you to find everything with extension in subdirs
file(GLOB sources *.cpp)

add_executable(express ${sources})

target_link_libraries(express ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} "${CMAKE_CURRENT_SOURCE_DIR}/../bamtools/lib/libbamtools.a" "pthread")

install(TARGETS express DESTINATION bin)
