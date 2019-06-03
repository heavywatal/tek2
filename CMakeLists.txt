cmake_minimum_required(VERSION 3.13)
project(teaposon
  VERSION 0.1.0
  LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)
include(GNUInstallDirs)
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
if(${CMAKE_INSTALL_PREFIX} MATCHES linuxbrew)
  set(CMAKE_INSTALL_LIBDIR lib)
endif()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_FLAGS_DEV "-O2 -g")
add_compile_options(-march=native -Wall -Wextra -pedantic)

set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

find_package(Threads)
find_package(ZLIB)
find_package(sfmt 0.7.2)
find_package(clippson 0.8.1)
find_package(wtl 0.8.4)

add_subdirectory(src)

add_executable(${PROJECT_NAME}-exe src/main.cpp)
target_link_libraries(${PROJECT_NAME}-exe PRIVATE objlib)
set_target_properties(${PROJECT_NAME}-exe PROPERTIES
  OUTPUT_NAME ${PROJECT_NAME}
)
install(TARGETS ${PROJECT_NAME}-exe
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()
