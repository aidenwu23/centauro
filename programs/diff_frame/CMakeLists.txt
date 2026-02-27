add_executable(df_match_jets match_jets.cc)
target_compile_features(df_match_jets PRIVATE cxx_std_17)
target_sources(
  df_match_jets
  PRIVATE
    ${CMAKE_SOURCE_DIR}/jet_tools/src/beam_helpers.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/electron_veto.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/math_helpers.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/transform_breit.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/root_io.cc)
target_include_directories(df_match_jets PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(
  df_match_jets
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4EIC::edm4eic
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::MathCore
    ROOT::GenVector
    ROOT::Tree
    ROOT::Hist
    ROOT::Gpad
    ROOT::Graf)

add_executable(splitting splitting.cc)
target_compile_features(splitting PRIVATE cxx_std_17)
target_sources(
  splitting
  PRIVATE
    ${CMAKE_SOURCE_DIR}/jet_tools/src/root_io.cc)
target_include_directories(splitting PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(
  splitting
  PRIVATE
    ROOT::Core
    ROOT::RIO
    ROOT::MathCore
    ROOT::GenVector
    ROOT::Tree
    ROOT::Hist
    ROOT::Gpad
    ROOT::Graf)

add_executable(eta_phi eta_phi.cc)
target_compile_features(eta_phi PRIVATE cxx_std_17)
target_sources(
  eta_phi
  PRIVATE
    ${CMAKE_SOURCE_DIR}/jet_tools/src/root_io.cc)
target_include_directories(eta_phi PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(
  eta_phi
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4EIC::edm4eic
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
    df_match_jets
    splitting
    eta_phi
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
