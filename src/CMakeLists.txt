# GLOB_RECURSE will allow you to find everything with extension in subdirs
file(GLOB sources *.cpp)

add_executable(cuffexpress ${sources})

target_link_libraries(cuffexpress ${Boost_LIBRARIES})

install(TARGETS cuffexpress DESTINATION bin)
