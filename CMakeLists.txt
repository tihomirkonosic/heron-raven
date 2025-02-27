cmake_minimum_required(VERSION 3.11)

project(raven VERSION 1.5.0
              LANGUAGES CXX
              DESCRIPTION "Raven is a de novo genome assembler for long uncorrected reads.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(FetchContent)
include(GNUInstallDirs)

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(raven_main_project ON)
endif ()
option(raven_build_tests "Build unit tests" ${raven_main_project})

find_package(ZLIB 1.2.8 REQUIRED)

find_package(cereal 1.3.0 QUIET)
if (NOT cereal_FOUND)
  FetchContent_Declare(
    cereal
    GIT_REPOSITORY https://github.com/USCiLab/cereal
    GIT_TAG v1.3.0)

  FetchContent_GetProperties(cereal)
  if (NOT cereal_POPULATED)
    FetchContent_Populate(cereal)
    add_subdirectory(
      ${cereal_SOURCE_DIR}
      ${cereal_BINARY_DIR}
      EXCLUDE_FROM_ALL)
    add_library(cereal::cereal ALIAS cereal)
  endif ()
endif ()

find_package(edlib 1.2.7 QUIET)
if (NOT edlib_FOUND)
  FetchContent_Declare(
    edlib
    GIT_REPOSITORY https://github.com/martinsos/edlib
    GIT_TAG v1.2.7)

  FetchContent_GetProperties(edlib)
  if (NOT edlib_POPULATED)
    FetchContent_Populate(edlib)
    add_subdirectory(
      ${edlib_SOURCE_DIR}
      ${edlib_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

find_package(thread_pool 4.0.0 QUIET)
if (NOT thread_pool_FOUND)
  FetchContent_Declare(
    thread_pool
    GIT_REPOSITORY https://github.com/rvaser/thread_pool
    GIT_TAG 4.0.0)

  FetchContent_GetProperties(thread_pool)
  if (NOT thread_pool_POPULATED)
    FetchContent_Populate(thread_pool)
    add_subdirectory(
      ${thread_pool_SOURCE_DIR}
      ${thread_pool_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

if (raven_build_tests)
  find_package(GTest 1.14.0 QUIET)
  if (NOT GTest_FOUND)
    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest
      GIT_TAG v1.14.0)

    FetchContent_GetProperties(googletest)
    if (NOT googletest_POPULATED)
      FetchContent_Populate(googletest)
      add_subdirectory(
        ${googletest_SOURCE_DIR}
        ${googletest_BINARY_DIR}
        EXCLUDE_FROM_ALL)
      add_library(GTest::Main ALIAS gtest_main)
    endif ()
  endif ()
endif ()

add_library(raven
  src/graph.cpp
  src/pile.cpp
  src/node.cpp
  src/edge.cpp
  src/parser.cpp
  src/graph_constructor/graph_constructor.cpp
  src/graph_assembler/graph_assembler.cpp
  src/option_manager.cpp
  src/graph_constructor/overlap_parser.cpp
  src/ram/minimizer_engine.cpp
  #src/extended_overlap.cpp
  #src/read_preprocessing.cpp
)

target_link_libraries(raven
  cereal::cereal
  edlib::edlib
  thread_pool::thread_pool
  ZLIB::ZLIB)

target_include_directories(raven PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

add_executable(raven_exe
  src/main.cpp
)

target_link_libraries(raven_exe
  raven)


target_compile_definitions(raven PUBLIC RAVEN_VERSION="${PROJECT_VERSION}")
set_property(TARGET raven_exe PROPERTY OUTPUT_NAME raven)

install(TARGETS raven_exe DESTINATION ${CMAKE_INSTALL_BINDIR})

if (raven_build_tests)
  enable_testing()

  add_executable(raven_test
    test/raven_test.cpp)

  target_link_libraries(raven_test
    raven
    GTest::Main)

  target_compile_definitions(raven_test
    PRIVATE TEST_DATA="${PROJECT_SOURCE_DIR}/test/data/")

  gtest_add_tests(TARGET raven_test SOURCES test/raven_test.cpp)
endif ()
