// Merge EDM4eic ROOT files by appending all events into one output.
/**
Run:
./run.sh programs/data_prep/merge_data.cc -o data/events/filtered/merged.root -i data/events/filtered/0000_skim.root,data/events/filtered/0001_skim.root
./run.sh programs/data_prep/merge_data.cc -o data/events/filtered/merged.root -i "data/events/filtered/*_skim.root" --drop-inputs
  
**/


#include <podio/Frame.h>
#include <podio/ROOTReader.h>
#include <podio/ROOTWriter.h>

#include <cstdlib>
#include <cctype>
#include <filesystem>
#include <glob.h>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace {

// CLI inputs: target output, list of inputs, optional cleanup.
struct Args {
  std::string output;
  std::vector<std::string> inputs;
  bool drop_inputs = false;
};

// CLI parser: expands comma-delimited inputs and optional drop flag.
void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -o OUTPUT.root -i INPUT.root [-i INPUT2.root ...]\n";
  std::cerr << "       Inputs may be glob patterns (quoted) or comma separated per -i.\n";
  std::cerr << "       Optional: --drop-inputs deletes merged inputs.\n";
}

Args parse_args(int argc, char* argv[]) {
  auto trim = [](std::string_view v) {
    while (!v.empty() && std::isspace(static_cast<unsigned char>(v.front()))) v.remove_prefix(1);
    while (!v.empty() && std::isspace(static_cast<unsigned char>(v.back()))) v.remove_suffix(1);
    return v;
  };

  Args args;
  for (int i = 1; i < argc; ++i) {
    std::string_view v(argv[i]);
    if ((v == "-o" || v == "--output") && i + 1 < argc) {
      args.output = argv[++i];
    } else if (v == "-i" || v == "--input") {
      if (i + 1 >= argc) {
        usage(argv[0]);
        std::exit(1);
      }
      std::string_view raw(argv[++i]);
      std::size_t start = 0;
      while (start < raw.size()) {
        const auto end = raw.find(',', start);
        std::string_view token = raw.substr(start, end == std::string_view::npos ? raw.size() - start : end - start);
        token = trim(token);
        if (!token.empty()) args.inputs.emplace_back(token);
        if (end == std::string_view::npos) break;
        start = end + 1;
      }
    } else if (v == "--drop-inputs") {
      args.drop_inputs = true;
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }
  if (args.output.empty() || args.inputs.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  return args;
}

}  // namespace

int main(int argc, char* argv[]) {
  const Args args = parse_args(argc, argv);

  // Expand glob patterns into explicit input lists before merging.
  std::vector<std::string> expanded;
  for (const auto& pat : args.inputs) {
    glob_t g{};
    const int res = glob(pat.c_str(), 0, nullptr, &g);
    if (res == GLOB_NOMATCH) {
      std::cerr << "[merge_data] No matches for pattern " << pat << "; skipping\n";
      globfree(&g);
      continue;
    }
    if (res != 0) {
      std::cerr << "[merge_data] glob failed for pattern " << pat << "\n";
      globfree(&g);
      continue;
    }
    for (std::size_t i = 0; i < g.gl_pathc; ++i) {
      expanded.emplace_back(g.gl_pathv[i]);
    }
    globfree(&g);
  }
  if (expanded.empty()) {
    std::cerr << "[merge_data] No inputs after glob expansion.\n";
    return 1;
  }

  podio::ROOTWriter writer(args.output);
  std::size_t files_merged = 0;
  std::size_t events_written = 0;

  // Append frames from each input file into the output writer.
  for (const auto& path : expanded) {
    podio::ROOTReader reader;
    try {
      reader.openFile(path);
    } catch (const std::exception& ex) {
      std::cerr << "[merge_data] Failed to open input " << path << ": " << ex.what() << "\n";
      continue;
    }

    const auto entries = reader.getEntries("events");
    if (entries == 0) {
      std::cerr << "[merge_data] No events in " << path << "; skipping\n";
      continue;
    }

    for (std::size_t i = 0; i < entries; ++i) {
      auto data = reader.readEntry("events", i);
      if (!data) continue;
      podio::Frame frame(std::move(data));
      writer.writeFrame(frame, "events");
      ++events_written;
    }
    ++files_merged;
  }

  writer.finish();

  std::cout << "[merge_data] Merged " << files_merged << " files into " << args.output << "\n";
  std::cout << "[merge_data] Wrote " << events_written << " events\n";

  // Optional cleanup: drop inputs after a successful merge.
  if (args.drop_inputs && !expanded.empty()) {
    const auto out_path = std::filesystem::absolute(args.output);
    for (const auto& in : expanded) {
      const auto in_path = std::filesystem::absolute(in);
      if (in_path == out_path) continue;
      std::error_code ec;
      std::filesystem::remove(in_path, ec);
      if (ec) {
        std::cerr << "[merge_data] Warning: failed to delete " << in << "\n";
      }
    }
  }
  return 0;
}
