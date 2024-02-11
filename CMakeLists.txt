cmake_minimum_required(VERSION 3.15)
project(tek2
  VERSION 0.2.1
  LANGUAGES CXX)

cmake_policy(SET CMP0076 NEW)
set(CMAKE_VERBOSE_MAKEFILE ON)
include(CMakePrintHelpers)
include(FetchContent)
include(GNUInstallDirs)
cmake_print_variables(CMAKE_INSTALL_PREFIX)
if(${CMAKE_INSTALL_PREFIX} MATCHES linuxbrew)
  set(CMAKE_INSTALL_LIBDIR lib)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
cmake_print_variables(CMAKE_BUILD_TYPE)
set(CMAKE_CXX_FLAGS_DEV "-O2 -g")
add_compile_options(-march=native -Wall -Wextra -pedantic)

set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

function(find_package_or_fetch package version repo)
  find_package(${package} ${version} QUIET)  # 3.24
  if(${package}_FOUND)
    cmake_print_variables(${package}_DIR)
    return()
  endif()
  string(TOUPPER "${package}" upackage)
  import_env("FETCHCONTENT_SOURCE_DIR_${upackage}")
  if(IS_DIRECTORY "${FETCHCONTENT_SOURCE_DIR_${upackage}}")
    set(tag "HEAD")
  else()
    set(tag "v${version}")
  endif()
  cmake_print_variables(FETCHCONTENT_SOURCE_DIR_${upackage} tag)
  set(FETCHCONTENT_QUIET OFF)
  FetchContent_Declare(
    ${package}
    GIT_REPOSITORY https://github.com/${repo}.git
    GIT_TAG ${tag}
    GIT_SHALLOW ON
    # FIND_PACKAGE_ARGS  # 3.24
  )
  # FetchContent_MakeAvailable(${package})  # 3.28
  if(NOT ${package}_POPULATED)
    FetchContent_Populate(${package})
    add_subdirectory(${${package}_SOURCE_DIR} ${${package}_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
  cmake_print_variables(${package}_SOURCE_DIR ${package}_BINARY_DIR)
endfunction()

find_package(Threads)
find_package(ZLIB)
find_package_or_fetch(wtl 0.9.0 heavywatal/cxxwtl)
find_package_or_fetch(clippson 0.8.6 heavywatal/clippson)
find_package_or_fetch(sfmt 0.7.3 heavywatal/sfmt-class)

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
