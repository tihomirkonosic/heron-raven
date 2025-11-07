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
include(ExternalProject)

# ---------------------------------------------------------
# Options
# ---------------------------------------------------------
if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(raven_main_project ON)
endif ()
option(raven_build_tests "Build unit tests" ${raven_main_project})

option(RAVEN_WITH_CATBOOST "Enable CatBoost integration" OFF)
set(RAVEN_CATBOOST_MODE "auto" CACHE STRING "catboost mode: auto|evaluator|header")
set_property(CACHE RAVEN_CATBOOST_MODE PROPERTY STRINGS auto evaluator header)

# Helper: find conan in user env (no admin needed)
function(_raven_find_conan outvar)
  find_program(_CONAN_EXECUTABLE NAMES conan conan.exe)
  if(_CONAN_EXECUTABLE)
    set(${outvar} "${_CONAN_EXECUTABLE}" PARENT_SCOPE)
  else()
    set(${outvar} "" PARENT_SCOPE)
  endif()
endfunction()

# ---------------------------------------------------------
# Dependencies
# ---------------------------------------------------------
find_package(biosoup 0.10.0 QUIET)
if (NOT biosoup_FOUND)
  FetchContent_Declare(
    biosoup
    GIT_REPOSITORY https://github.com/rvaser/biosoup
    GIT_TAG 0.10.0)
  FetchContent_GetProperties(biosoup)
  if (NOT biosoup_POPULATED)
    FetchContent_Populate(biosoup)
    add_subdirectory(${biosoup_SOURCE_DIR} ${biosoup_BINARY_DIR} EXCLUDE_FROM_ALL)
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
    add_subdirectory(${bioparser_SOURCE_DIR} ${bioparser_BINARY_DIR} EXCLUDE_FROM_ALL)
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
    add_subdirectory(${cereal_SOURCE_DIR} ${cereal_BINARY_DIR} EXCLUDE_FROM_ALL)
    add_library(cereal::cereal ALIAS cereal)
  endif ()
endif ()

find_package(ram 1.0.0 QUIET)
if (NOT ram_FOUND)
  FetchContent_Declare(
    ram
    GIT_REPOSITORY https://github.com/FilipTomas/ram_count
    GIT_TAG master
    GIT_SUBMODULES "")
  FetchContent_GetProperties(ram)
  if (NOT ram_POPULATED)
    FetchContent_Populate(ram)
    add_subdirectory(${ram_SOURCE_DIR} ${ram_BINARY_DIR} EXCLUDE_FROM_ALL)
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
    add_subdirectory(${edlib_SOURCE_DIR} ${edlib_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif ()
endif ()

# ---------------------------------------------------------
# CatBoost integration (evaluator or header-only fallback)
# ---------------------------------------------------------
set(RAVEN_CATBOOST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/third_party/catboost_model/catboost_model.hpp")
set(RAVEN_CATBOOST_TARGET "")

if (RAVEN_WITH_CATBOOST)
  set(_USE_EVALUATOR FALSE)

  if (RAVEN_CATBOOST_MODE STREQUAL "evaluator" OR RAVEN_CATBOOST_MODE STREQUAL "auto")
    _raven_find_conan(_CONAN)
    if (_CONAN)
      FetchContent_Declare(
        catboost_src
        GIT_REPOSITORY https://github.com/catboost/catboost.git
        GIT_TAG v1.2.5
        GIT_SHALLOW TRUE)
      FetchContent_GetProperties(catboost_src)
      if (NOT catboost_src_POPULATED)
        FetchContent_Populate(catboost_src)
      endif ()

      set(CATBOOST_BUILD_DIR   ${CMAKE_BINARY_DIR}/_dep_catboost/build)
      set(CATBOOST_INSTALL_DIR ${CMAKE_BINARY_DIR}/_dep_catboost/install)

      ExternalProject_Add(catboost_build
        SOURCE_DIR  ${catboost_src_SOURCE_DIR}
        BINARY_DIR  ${CATBOOST_BUILD_DIR}
        INSTALL_DIR ${CATBOOST_INSTALL_DIR}
        CMAKE_ARGS
          -DCMAKE_BUILD_TYPE=$<IF:$<CONFIG:>,Release,$<CONFIG>>
          -DCMAKE_CXX_STANDARD=17
          -DCMAKE_INSTALL_PREFIX=${CATBOOST_INSTALL_DIR}
          -DBUILD_JAVA=OFF -DBUILD_PYTHON=OFF -DBUILD_R=OFF
          -DUSE_CONAN=ON
          -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
          -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        BUILD_BYPRODUCTS
          ${CATBOOST_BUILD_DIR}/lib/libcatboostmodel.a
      )

      add_library(catboostmodel STATIC IMPORTED GLOBAL)
      add_dependencies(catboostmodel catboost_build)

      if (APPLE)
        set(_CATBOOST_LIB ${CATBOOST_BUILD_DIR}/lib/libcatboostmodel.a)
      elseif (WIN32)
        set(_CATBOOST_LIB ${CATBOOST_BUILD_DIR}/lib/catboostmodel.lib)
      else ()
        set(_CATBOOST_LIB ${CATBOOST_BUILD_DIR}/lib/libcatboostmodel.a)
      endif ()

      set_target_properties(catboostmodel PROPERTIES
        IMPORTED_LOCATION ${_CATBOOST_LIB}
        INTERFACE_INCLUDE_DIRECTORIES
          "${catboost_src_SOURCE_DIR};${catboost_src_SOURCE_DIR}/catboost"
      )

      if (UNIX AND NOT APPLE)
        target_link_libraries(catboostmodel INTERFACE pthread dl)
      endif ()

      target_compile_definitions(catboostmodel INTERFACE RAVEN_CATBOOST_EVALUATOR=1)
      add_library(catboost::model ALIAS catboostmodel)

      set(RAVEN_CATBOOST_TARGET catboostmodel)
      set(_USE_EVALUATOR TRUE)
      message(STATUS "CatBoost evaluator enabled (Conan found).")
    else()
      message(WARNING "Conan not found; falling back to CatBoost header-only mode.")
    endif()
  endif()

  if (RAVEN_CATBOOST_MODE STREQUAL "header" OR NOT _USE_EVALUATOR)
    if (NOT EXISTS "${RAVEN_CATBOOST_HEADER}")
      message(FATAL_ERROR "CatBoost header '${RAVEN_CATBOOST_HEADER}' not found.
Export from Python:  model.save_model('third_party/catboost_model/catboost_model.hpp', format='cpp')")
    endif()

    add_library(catboost_header_only INTERFACE)
    target_include_directories(catboost_header_only INTERFACE
      ${CMAKE_CURRENT_SOURCE_DIR}/third_party/catboost_model)
    target_compile_definitions(catboost_header_only INTERFACE RAVEN_CATBOOST_HEADER=1)
    add_library(catboost::model ALIAS catboost_header_only)

    set(RAVEN_CATBOOST_TARGET catboost_header_only)
    message(STATUS "CatBoost header-only mode active.")
  endif()
endif ()

# ---------------------------------------------------------
# Your targets
# ---------------------------------------------------------
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
  src/graph_constructor/ram_overlap.cpp
  #src/extended_overlap.cpp
  #src/read_preprocessing.cpp
)

# Use KEYWORD signature everywhere for this target
target_link_libraries(raven PUBLIC
  cereal::cereal
  biosoup::biosoup
  bioparser::bioparser
  ram::ram
  edlib::edlib
)

# Link CatBoost only if enabled
if (RAVEN_WITH_CATBOOST AND RAVEN_CATBOOST_TARGET)
  target_link_libraries(raven PUBLIC catboost::model)
endif ()

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

# ---------------------------------------------------------
# Tests
# ---------------------------------------------------------
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
      add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL)
      add_library(GTest::Main ALIAS gtest_main)
    endif ()
  endif ()

  enable_testing()
  add_executable(raven_test test/raven_test.cpp)
  target_link_libraries(raven_test raven bioparser::bioparser GTest::Main)
  target_compile_definitions(raven_test PRIVATE TEST_DATA="${PROJECT_SOURCE_DIR}/test/data/")
  gtest_add_tests(TARGET raven_test SOURCES test/raven_test.cpp)
endif ()
