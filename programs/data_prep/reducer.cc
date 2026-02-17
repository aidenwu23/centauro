// Reduce an EDM4eic file to a subset of collections.
/*
Run:
./run.sh programs/data_prep/reducer.cc -i data/events/raw/events.root -o data/events/filtered/events_skim.root
./run.sh programs/data_prep/reducer.cc -i data/events/raw/events.root --list
*/

#include <podio/Frame.h>
#include <podio/ROOTReader.h>
#include <podio/ROOTWriter.h>

#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace {

// Collections kept in the skim output.
const std::vector<std::string> kKeep = {
    "EventHeader",
    "MCParticles",
    "GeneratedParticles",
    "GeneratedBreitFrameParticles",
    "ReconstructedParticles",
    "ReconstructedBreitFrameParticles",
    // "GeneratedJets",
    // "GeneratedCentauroJets",
    // "ReconstructedJets",
    // "ReconstructedCentauroJets",
    "ReconstructedChargedParticles",
    "ReconstructedChargedParticleAssociations",
    "ReconstructedElectrons",
    "ReconstructedElectronsForDIS",
    "ScatteredElectronsTruth",
    "InclusiveKinematicsTruth",
    "InclusiveKinematicsElectron",
    // "InclusiveKinematicsJB",
    // "InclusiveKinematicsDA",
    "ReconstructedParticleAssociations",
    // "HadronicFinalState",
    // Track collections needed for ReconstructedParticle.getTracks() and measurements_size().
    "CentralCKFTracks",
    "CentralCKFTracksUnfiltered",
    "B0TrackerCKFTracks",
    "B0TrackerCKFTracksUnfiltered",
    "CentralTrackerMeasurements",
    "B0TrackerMeasurements",
};

// Parsed command-line inputs.
struct Args {
  std::string input;
  std::string output;
  bool list_only = false;
};

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT.root [-o OUTPUT.root] [--list]\n";
}

// Parse minimal CLI flags for this tool.
Args parse_args(int argc, char* argv[]) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    std::string_view v(argv[i]);
    if ((v == "-i" || v == "--input") && i + 1 < argc) {
      args.input = argv[++i];
    } else if ((v == "-o" || v == "--output") && i + 1 < argc) {
      args.output = argv[++i];
    } else if (v == "--list") {
      args.list_only = true;
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }
  if (args.input.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  if (!args.list_only && args.output.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  return args;
}

}  // namespace

int main(int argc, char* argv[]) {
  const Args args = parse_args(argc, argv);

  // Open the input file and fail fast on errors.
  podio::ROOTReader reader;
  try {
    reader.openFile(args.input);
  } catch (const std::exception& ex) {
    std::cerr << "[reducer] Failed to open input: " << ex.what() << "\n";
    return 1;
  }

  const auto entries = reader.getEntries("events");
  // Stop early when the input carries no events.
  if (entries == 0) {
    std::cerr << "[reducer] No events in input.\n";
    return 1;
  }

  // List available collections if requested.
  if (args.list_only) {
    auto data = reader.readEntry("events", 0);
    if (!data) {
      std::cerr << "[reducer] Failed to read first event for listing.\n";
      return 1;
    }
    podio::Frame frame(std::move(data));
    auto names = frame.getAvailableCollections();
    std::cout << "[reducer] Collections in input:\n";
    for (const auto& name : names) {
      std::cout << "  " << name << "\n";
    }
    return 0;
  }

  podio::ROOTWriter writer(args.output);

  std::size_t written = 0;
  std::size_t missing_warns = 0;
  // Copy only the collections we keep for each event.
  for (std::size_t i = 0; i < entries; ++i) {
    auto data = reader.readEntry("events", i);
    if (!data) continue;
    podio::Frame frame(std::move(data));
    // Collect available names once to make membership checks fast.
    std::unordered_set<std::string> present_names;
    for (const auto& name : frame.getAvailableCollections()) {
      present_names.emplace(name);
    }

    std::vector<std::string> present;
    present.reserve(kKeep.size());
    for (const auto& name : kKeep) {
      if (present_names.count(name)) present.push_back(name);
    }
    if (present.empty()) {
      if (missing_warns < 5) {
        std::cerr << "[reducer] Warning: no requested collections present in event " << i << "\n";
      }
      ++missing_warns;
      continue;
    }

    writer.writeFrame(frame, "events", present);
    ++written;
  }

  writer.finish();

  std::cout << "[reducer] Wrote " << written << " events to " << args.output << "\n";
  std::cout << "[reducer] Kept " << kKeep.size() << " collections.\n";
  if (missing_warns > 0) {
    std::cout << "[reducer] Events with no kept collections: " << missing_warns << "\n";
  }
  return 0;
}
