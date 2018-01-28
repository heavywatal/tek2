cmake_minimum_required(VERSION 3.1)
project(tek CXX)

include(GNUInstallDirs)
option(PROFILE "Use gperftools" FALSE)
set(CMAKE_VERBOSE_MAKEFILE TRUE)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_compile_options(-g)
else()
  add_definitions(-DNDEBUG)
endif()
add_compile_options(-O3 -march=native -Wall -Wextra -Wpedantic)

find_package(wtl)
find_package(sfmt)

set(Boost_NO_BOOST_CMAKE TRUE)
find_package(Boost REQUIRED COMPONENTS program_options iostreams filesystem system)
message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")

message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
if(${CMAKE_INSTALL_PREFIX} MATCHES linuxbrew)
  set(CMAKE_INSTALL_LIBDIR lib)
endif()

set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND git describe --always --abbrev=0 --dirty
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND git show -s --format=%ci
  OUTPUT_VARIABLE GIT_COMMIT_TIME
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "commit ${GIT_COMMIT_HASH} [${GIT_BRANCH}]")
message(STATUS "Date:  ${GIT_COMMIT_TIME}")
configure_file(
  ${CMAKE_SOURCE_DIR}/src/version.hpp.in
  ${CMAKE_SOURCE_DIR}/src/version.hpp @ONLY
)

aux_source_directory(src srcfiles)
message(STATUS "Source files: ${srcfiles}")

add_library(${PROJECT_NAME} SHARED ${srcfiles})
target_include_directories(${PROJECT_NAME} BEFORE
  PRIVATE ${Boost_INCLUDE_DIRS}
  INTERFACE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_link_libraries(${PROJECT_NAME} PRIVATE wtl sfmt ${Boost_LIBRARIES})

add_executable(${PROJECT_NAME}-exe main.cpp)
target_link_libraries(${PROJECT_NAME}-exe ${PROJECT_NAME})
if(PROFILE)
  target_link_libraries(${PROJECT_NAME}-exe profiler)
endif()
set_target_properties(${PROJECT_NAME}-exe PROPERTIES
  OUTPUT_NAME ${PROJECT_NAME}
)

install(TARGETS ${PROJECT_NAME}-exe ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}-config
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(DIRECTORY src/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  FILES_MATCHING PATTERN "*.hpp")
install(EXPORT ${PROJECT_NAME}-config
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}
)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()
