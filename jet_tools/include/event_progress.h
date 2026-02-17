#pragma once

#include <RtypesCore.h>

#include <cstddef>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_set>

namespace jet_tools {

class EventLimitGate {
 public:
  explicit EventLimitGate(std::size_t max_events)
      : max_events_(max_events),
        limited_(max_events != std::numeric_limits<std::size_t>::max()) {}

  bool allow(ULong64_t event_id) {
    // Accept all rows belonging to the first max_events unique event IDs in stream order.
    if (!limited_) {
      return true;
    }
    if (max_events_ == 0) {
      return false;
    }
    if (accepted_event_ids_.find(event_id) != accepted_event_ids_.end()) {
      return true;
    }
    if (accepted_event_ids_.size() >= max_events_) {
      return false;
    }
    accepted_event_ids_.insert(event_id);
    return true;
  }

 private:
  std::size_t max_events_ = std::numeric_limits<std::size_t>::max();
  bool limited_ = false;
  std::unordered_set<ULong64_t> accepted_event_ids_;
};

class ProgressTicker {
 public:
  explicit ProgressTicker(std::size_t report_every = 100)
      : report_every_(report_every) {}

  bool should_report(std::size_t index) const {
    return report_every_ > 0 && index % report_every_ == 0;
  }

  void report(const std::string& message) {
    // Overwrite the current terminal line and clear leftover characters from older text.
    std::cout << '\r' << message;
    if (last_chars_ > message.size()) {
      std::cout << std::string(last_chars_ - message.size(), ' ');
    }
    std::cout << std::flush;
    last_chars_ = message.size();
  }

  void finish(const std::string& message) {
    report(message);
    std::cout << "\n";
  }

 private:
  std::size_t report_every_ = 100;
  std::size_t last_chars_ = 0;
};

class EventCounter {
 public:
  bool observe(ULong64_t event_id) {
    // Count a new event only when event_id changes while streaming rows.
    if (have_last_ && event_id == last_event_id_) {
      return false;
    }
    last_event_id_ = event_id;
    have_last_ = true;
    ++events_seen_;
    return true;
  }

  std::size_t events_seen() const {
    return events_seen_;
  }

 private:
  ULong64_t last_event_id_ = 0;
  bool have_last_ = false;
  std::size_t events_seen_ = 0;
};

}  // namespace jet_tools
