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

find_package(biosoup 0.10.0 QUIET)
if (NOT biosoup_FOUND)
  FetchContent_Declare(
    biosoup
    GIT_REPOSITORY https://github.com/rvaser/biosoup
    GIT_TAG 0.10.0)

  FetchContent_GetProperties(biosoup)
  if (NOT biosoup_POPULATED)
    FetchContent_Populate(biosoup)
    add_subdirectory(
      ${biosoup_SOURCE_DIR}
      ${biosoup_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

find_package(bioparser 3.0.13 QUIET)
if (NOT bioparser_FOUND)
  FetchContent_Declare(
    bioparser
    GIT_REPOSITORY https://github.com/rvaser/bioparser
    GIT_TAG 3.0.13)

  FetchContent_GetProperties(bioparser)
  if (NOT bioparser_POPULATED)
    FetchContent_Populate(bioparser)
    add_subdirectory(
      ${bioparser_SOURCE_DIR}
      ${bioparser_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

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

find_package(ram 1.0.0 QUIET)
if (NOT ram_FOUND)
  FetchContent_Declare(
    ram
    GIT_REPOSITORY https://github.com/kocijan/ram
    GIT_TAG kocijan/hpc)

  FetchContent_GetProperties(ram)
  if (NOT ram_POPULATED)
    FetchContent_Populate(ram)
    add_subdirectory(
      ${ram_SOURCE_DIR}
      ${ram_BINARY_DIR}
      EXCLUDE_FROM_ALL)
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
  src/overlap_parser.cpp
  src/polymer_manager.cpp
  src/graph_assembler/graph_assembler_transitive.cpp
  src/graph_assembler/graph_assembler_tips.cpp
  src/graph_assembler/graph_assembler_bubbles.cpp
  #src/extended_overlap.cpp
  #src/read_preprocessing.cpp
)

target_link_libraries(raven
  cereal::cereal
  biosoup::biosoup
  bioparser::bioparser
  ram::ram
  edlib::edlib)

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
    bioparser::bioparser
    GTest::Main)

  target_compile_definitions(raven_test
    PRIVATE TEST_DATA="${PROJECT_SOURCE_DIR}/test/data/")

  gtest_add_tests(TARGET raven_test SOURCES test/raven_test.cpp)
endif ()
