add_executable(reducer reducer.cc)
target_compile_features(reducer PRIVATE cxx_std_17)
target_include_directories(reducer PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(
  reducer
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    EDM4EIC::edm4eic
    ROOT::Core
    ROOT::RIO)

add_executable(cut_events cut_events.cc)
target_compile_features(cut_events PRIVATE cxx_std_17)
target_sources(
  cut_events
  PRIVATE
    ${CMAKE_SOURCE_DIR}/jet_tools/src/beam_helpers.cc
    ${CMAKE_SOURCE_DIR}/jet_tools/src/math_helpers.cc)
target_include_directories(
  cut_events
  PRIVATE
    ${CMAKE_SOURCE_DIR})
target_link_libraries(
  cut_events
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    EDM4EIC::edm4eic
    ROOT::Core
    ROOT::RIO
    ROOT::MathCore
    ROOT::GenVector)

add_executable(merge_data merge_data.cc)
target_compile_features(merge_data PRIVATE cxx_std_17)
target_include_directories(merge_data PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(
  merge_data
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    EDM4EIC::edm4eic
    ROOT::Core
    ROOT::RIO)

install(
  TARGETS
    reducer
    cut_events
    merge_data
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
