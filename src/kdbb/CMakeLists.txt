add_subdirectory(pmc)

add_library(KDBB kdbb.h kdbb.cpp)
target_link_libraries(KDBB Utils PMC)
