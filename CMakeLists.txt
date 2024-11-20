cmake_minimum_required(VERSION 3.15)
project(tek2
  VERSION 0.2.2
  LANGUAGES CXX)

cmake_policy(SET CMP0169 OLD)  # 3.28
set(CMAKE_VERBOSE_MAKEFILE ON)
include(CMakePrintHelpers)
include(FetchContent)
include(GNUInstallDirs)
cmake_print_variables(CMAKE_INSTALL_PREFIX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
cmake_print_variables(CMAKE_BUILD_TYPE)
set(CMAKE_CXX_FLAGS_DEV "-O2 -g")

function(import_env variable)
  if(DEFINED ENV{${variable}})
    set(${variable} $ENV{${variable}} PARENT_SCOPE)
  endif()
endfunction()

function(find_package_or_fetch package tag repo)
  string(REGEX REPLACE "^v" "" version ${tag})
  find_package(${package} ${version} QUIET)  # 3.24
  if(${package}_FOUND)
    cmake_print_variables(${package}_DIR)
    return()
  endif()
  string(TOUPPER "${package}" upackage)
  import_env("FETCHCONTENT_SOURCE_DIR_${upackage}")
  if(IS_DIRECTORY "${FETCHCONTENT_SOURCE_DIR_${upackage}}")
    set(tag "HEAD")
  endif()
  cmake_print_variables(FETCHCONTENT_SOURCE_DIR_${upackage} tag)
  set(FETCHCONTENT_QUIET OFF)
  FetchContent_GetProperties(${package})
  FetchContent_Declare(
    ${package}
    GIT_REPOSITORY https://github.com/${repo}.git
    GIT_TAG ${tag}
    GIT_SHALLOW ON
    # FIND_PACKAGE_ARGS  # 3.24
  )
  # FetchContent_MakeAvailable(${package})  # 3.28
  FetchContent_GetProperties(${package})
  if(NOT ${package}_POPULATED)
    FetchContent_Populate(${package})
    add_subdirectory(${${package}_SOURCE_DIR} ${${package}_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
  cmake_print_variables(${package}_SOURCE_DIR ${package}_BINARY_DIR)
endfunction()

find_package(Threads)
find_package_or_fetch(wtl v0.9.2 heavywatal/cxxwtl)
find_package_or_fetch(clippson v0.8.8 heavywatal/clippson)
find_package_or_fetch(pcglite v0.1.2 heavywatal/pcglite)

add_library(${PROJECT_NAME} STATIC)
target_compile_options(${PROJECT_NAME} PUBLIC
  -Wall -Wextra -pedantic
  $<$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},x86_64>:-march=native>
  $<$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},arm64>:-march=armv8.3-a+sha3>
)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
)
target_link_libraries(${PROJECT_NAME} PUBLIC
  wtl::wtl
  wtl::zlib
  pcglite::pcglite
  Threads::Threads
  PRIVATE clippson::clippson
)
add_subdirectory(src)

add_executable(${PROJECT_NAME}-exe src/main.cpp)
target_link_libraries(${PROJECT_NAME}-exe PRIVATE ${PROJECT_NAME})
set_target_properties(${PROJECT_NAME}-exe PROPERTIES
  OUTPUT_NAME ${PROJECT_NAME}
)
install(TARGETS ${PROJECT_NAME}-exe)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()
