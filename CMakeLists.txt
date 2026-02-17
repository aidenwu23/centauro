cmake_minimum_required(VERSION 3.24)
cmake_policy(SET CMP0144 NEW)

project(centauro VERSION 0.1.0 LANGUAGES CXX)

# Build settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Keep build outputs in the build tree; installs still go under CMAKE_INSTALL_PREFIX.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

set(default_build_type "RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}'")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(DEFINED ENV{EIC_SHELL_PREFIX})
  list(PREPEND CMAKE_PREFIX_PATH "$ENV{EIC_SHELL_PREFIX}")
  list(APPEND CMAKE_MODULE_PATH "$ENV{EIC_SHELL_PREFIX}/lib/cmake/EICrecon")
endif()
list(APPEND CMAKE_PREFIX_PATH "/opt/local")
# Pick up common Geant4 modules for optional packages.
list(APPEND CMAKE_MODULE_PATH "/opt/local/lib/cmake/Geant4/Modules")

include(GNUInstallDirs)

# RPATH settings similar to EICrecon
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
set(CMAKE_SKIP_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

option(CENTAURO_WITH_EDM "Build EDM4hep/EDM4EIC tools" ON)

set(EDM4EIC_VERSION_MIN 5.0)
set(EDM4HEP_VERSION_MIN 0.7.1)
set(podio_VERSION_MIN 1.3.0)
set(ROOT_VERSION_MIN 6.28)
set(HepMC3_VERSION_MIN 3.2)

if(CENTAURO_WITH_EDM)
  # EIC environment provides CMake config packages for these
  find_package(podio ${podio_VERSION_MIN} REQUIRED)
  find_package(EDM4HEP ${EDM4HEP_VERSION_MIN} REQUIRED)
  find_package(EDM4EIC ${EDM4EIC_VERSION_MIN} REQUIRED)
  find_package(FastJet REQUIRED)
  find_package(FastJetTools REQUIRED)
  find_package(ROOT ${ROOT_VERSION_MIN} REQUIRED COMPONENTS Core RIO MathCore GenVector Tree Hist Gpad Graf TreePlayer)

  # fjcontrib Centauro plugin lives in the contrib "fragile" lib under the EIC env.
  find_library(
    FJCONTRIB_FRAGILE_LIB
    NAMES fastjetcontribfragile
    HINTS
      $ENV{EIC_SHELL_PREFIX}/lib)
  if(NOT FJCONTRIB_FRAGILE_LIB)
    message(FATAL_ERROR "fastjetcontribfragile (Centauro plugin) not found; ensure fjcontrib is available in the environment")
  endif()

  add_subdirectory(programs)
endif()
