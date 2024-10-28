execute_process(
  COMMAND git describe --always --dirty
  OUTPUT_VARIABLE GIT_DESCRIBE
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(GIT_DESCRIBE)
  set(PROJECT_VERSION ${GIT_DESCRIBE})
endif()
cmake_print_variables(PROJECT_VERSION)

configure_file(
  version.cpp.in
  version.cpp @ONLY
)

add_library(objlib OBJECT
  haploid.cpp
  population.cpp
  program.cpp
  transposon.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
)
target_compile_options(objlib PUBLIC
  -Wall -Wextra -pedantic
  $<$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},x86_64>:-march=native>
  $<$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},arm64>:-march=armv8.3-a+sha3>
)
target_compile_features(objlib PUBLIC cxx_std_17)
set_target_properties(objlib PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(objlib PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
)
target_link_libraries(objlib PUBLIC
  wtl::wtl
  wtl::sfmt
  clippson::clippson
  ZLIB::ZLIB
  Threads::Threads
)
