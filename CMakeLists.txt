cmake_minimum_required(VERSION 3.8)
project(tek2
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

function(find_package_or_git package version repo)
  find_package(${package} ${version} QUIET)
  if(${package}_FOUND)
    return()
  endif()
  set(subdir ${PROJECT_BINARY_DIR}/${package}-v${version})
  if(NOT EXISTS ${subdir})
    find_package(Git)
    execute_process(COMMAND
      ${GIT_EXECUTABLE} clone --recursive --depth=1 --branch=v${version} https://github.com/${repo}.git ${subdir}
    )
  endif()
  message(STATUS "add_subdirectory(${subdir})")
  add_subdirectory(${subdir} EXCLUDE_FROM_ALL)
endfunction()

find_package(Threads)
find_package(ZLIB)
find_package_or_git(wtl 0.8.5 heavywatal/cxxwtl)
find_package_or_git(clippson 0.8.3 heavywatal/clippson)
find_package_or_git(sfmt 0.7.2 heavywatal/sfmt-class)

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
