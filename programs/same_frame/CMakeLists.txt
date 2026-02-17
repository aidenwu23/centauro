add_executable(sf_match_jets match_jets.cc)
target_compile_features(sf_match_jets PRIVATE cxx_std_17)
target_sources(
  sf_match_jets
  PRIVATE
    ${CMAKE_SOURCE_DIR}/jet_tools/src/math_helpers.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/root_io.cc)
target_include_directories(sf_match_jets PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(
  sf_match_jets
  PRIVATE
    ROOT::Core
    ROOT::RIO
    ROOT::MathCore
    ROOT::GenVector
    ROOT::Tree
    ROOT::Hist
    ROOT::Gpad
    ROOT::Graf)

add_executable(performance performance.cc)
target_compile_features(performance PRIVATE cxx_std_17)
target_sources(
  performance
  PRIVATE
    ${CMAKE_SOURCE_DIR}/jet_tools/src/math_helpers.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/plot_helpers.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/root_io.cc)
target_include_directories(performance PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(
  performance
  PRIVATE
    ROOT::Core
    ROOT::RIO
    ROOT::MathCore
    ROOT::GenVector
    ROOT::Tree
    ROOT::Hist
    ROOT::Gpad
    ROOT::Graf)

add_executable(controls controls.cc)
target_compile_features(controls PRIVATE cxx_std_17)
target_sources(
  controls
  PRIVATE
    ${CMAKE_SOURCE_DIR}/jet_tools/src/root_io.cc)
target_include_directories(controls PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(
  controls
  PRIVATE
    ROOT::Core
    ROOT::RIO
    ROOT::MathCore
    ROOT::GenVector
    ROOT::Tree
    ROOT::Hist
    ROOT::Gpad
    ROOT::Graf)

install(
  TARGETS
    sf_match_jets
    performance
    controls
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
