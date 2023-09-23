cmake_minimum_required(VERSION 3.0.0)

include_directories(utils defclique)

add_subdirectory(utils)
add_subdirectory(defclique)

set(EXECUTABLE_NAME run)

add_executable(${EXECUTABLE_NAME} main.cpp)
target_link_libraries(${EXECUTABLE_NAME} Utils DefClique)
