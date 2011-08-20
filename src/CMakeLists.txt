# GLOB_RECURSE will allow you to find everything with extension in subdirs
file(GLOB sources *.cpp)

add_executable(express ${sources})

target_link_libraries(express ${Boost_LIBRARIES} "bamtools")

install(TARGETS express DESTINATION bin)
