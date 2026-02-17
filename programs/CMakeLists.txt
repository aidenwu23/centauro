set(CENTAURO_PROGRAM_INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}
  ${FASTJET_INCLUDE_DIRS}
  ${FASTJET_TOOLS_INCLUDE_DIRS})

set(CENTAURO_ROOT_LIBS
  ROOT::Core
  ROOT::RIO
  ROOT::MathCore
  ROOT::GenVector
  ROOT::Tree
  ROOT::Hist
  ROOT::Gpad
  ROOT::Graf)

set(CENTAURO_EDM_LIBS
  podio::podio
  podio::podioRootIO
  EDM4EIC::edm4eic
  EDM4HEP::edm4hep)

set(CENTAURO_FASTJET_LIBS
  ${FASTJET_LIBRARIES}
  ${FASTJET_TOOLS_LIBRARIES}
  ${FASTJET_LIBRARY}
  ${FJTOOLS_LIBRARY}
  ${FJCONTRIB_FRAGILE_LIB})

add_executable(cluster_jets cluster_jets.cc)
target_compile_features(cluster_jets PRIVATE cxx_std_17)
target_sources(
  cluster_jets
  PRIVATE
    ${CMAKE_SOURCE_DIR}/jet_tools/src/electron_veto.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/jet_builders.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/math_helpers.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/beam_helpers.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/root_io.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/transform_breit.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/z_jet.cc)
target_include_directories(cluster_jets PRIVATE ${CENTAURO_PROGRAM_INCLUDE_DIRS})
target_link_libraries(
  cluster_jets
  PRIVATE
    ${CENTAURO_EDM_LIBS}
    ${CENTAURO_ROOT_LIBS}
    ${CENTAURO_FASTJET_LIBS})

add_executable(eta_vs_z eta_vs_z.cc)
target_compile_features(eta_vs_z PRIVATE cxx_std_17)
target_include_directories(eta_vs_z PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(
  eta_vs_z
  PRIVATE
    ${CENTAURO_ROOT_LIBS})

add_executable(overlay overlay.cc)
target_compile_features(overlay PRIVATE cxx_std_17)
target_sources(
  overlay
  PRIVATE
    ${CMAKE_SOURCE_DIR}/jet_tools/src/plot_helpers.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/math_helpers.cc)
target_include_directories(overlay PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(overlay PRIVATE ${CENTAURO_ROOT_LIBS})

add_subdirectory(data_prep)
add_subdirectory(same_frame)
add_subdirectory(diff_frame)

install(
  TARGETS
    cluster_jets
    eta_vs_z
    overlay
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
