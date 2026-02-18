# centauro

This directory contains the C++ analysis pipeline used to compare anti-k_t and Centauro jets on EDM4eic/EDM4hep event data. The workflow is:

(0. Prepare events.)
1. Cluster jets.
2. Match jets.
3. Produce analysis ROOT outputs and overlays.

## Directory layout

- `build.sh`: configure, build, install helper for CMake.
- `run.sh`: run a built executable by source path or binary path.
- `CMakeLists.txt`: top-level build configuration.
- `programs/`: executable entry points.
- `jet_tools/`: shared helpers used by `programs/`.
- `tools/read_outputs.py`: pull summary numbers from ROOT outputs, probably outdated.
- `data/`: generated inputs/outputs.
  - `data/events/raw`
  - `data/events/filtered`
  - `data/jets/raw`
  - `data/jets/matched`
  - `data/graphs`

## Build

Run:

```bash
./build.sh
```

## Run helper

Use:

```bash
./run.sh <path-to-source-or-binary> [program args...]
```

Examples:

```bash
./run.sh programs/cluster_jets.cc -h
./run.sh build/cluster_jets -h
```

`run.sh` resolves source paths whose executable names differ from file basenames, including:

- `programs/same_frame/match_jets.cc` -> `build/sf_match_jets`
- `programs/diff_frame/match_jets.cc` -> `build/df_match_jets`

## Workflow

### 1) Cluster jets

Native-threshold output (Lab particles and Breit particles receive the same numerical cuts, which are non-invariant):

```bash
./run.sh programs/cluster_jets.cc \
  -i data/events/filtered/merged.root \
  -o data/jets/raw/jets_native.root \
  --threshold-frame native
```

Lab-threshold output (Lab particles are used to build a set of skip indices, which is mapped to the Breit particles so Breit particles get lab frame cuts as well):

```bash
./run.sh programs/cluster_jets.cc \
  -i data/events/filtered/merged.root \
  -o data/jets/raw/jets_lab.root \
  --threshold-frame lab
```

Defaults:

- Jet clustering radius: `-R 0.6`
- Pseudorapidity:`--max-eta 4`
- Minimum jet transverse momenta: `--min-jet-pT 5`
- Frame in which the kinematic cuts above are applied: `--threshold-frame native`

### 2A) Different-frame analysis

Cross-frame matching (take a Centauro jet clustered in Breit frame, boost it back to lab frame, and match to the closets lab-clustered anti-kt jet):

```bash
./run.sh programs/diff_frame/match_jets.cc \
  -j data/jets/raw/jets_lab.root \
  -e data/events/filtered/merged.root \
  -o data/jets/matched/diff_frame_matches.root \
  --match-dR 0.3
```

Splitting and overlap metrics:

```bash
./run.sh programs/diff_frame/splitting.cc \
  -m data/jets/matched/diff_frame_matches.root \
  -j data/jets/raw/jets_lab.root \
  -o data/graphs/diff_frame_splitting.root \
  --max-dR 0.3
```

### 2A) Same-frame analysis

Take jets clustered on truth particles and match to jets clustered on reco particles of the same frame ranked by transverse-momenta similarity.

```bash
./run.sh programs/same_frame/match_jets.cc \
  -i data/jets/raw/jets_native.root \
  -o data/jets/matched/same_frame_matches.root
```

Performance plots:

```bash
./run.sh programs/same_frame/performance.cc \
  -m data/jets/matched/same_frame_matches.root \
  -j data/jets/raw/jets_native.root \
  -o data/graphs/performance.root \
  --abs-eta-max 3 --pT-min 5
```

Control plots:

```bash
./run.sh programs/same_frame/controls.cc \
  -m data/jets/matched/same_frame_matches.root \
  -j data/jets/raw/jets_native.root \
  -o data/graphs/controls.root \
  --abs-eta-max 3 --pT-min 5 --match-dR 0.2
```

### 4) Additional plots

Eta vs z maps (see Centauro paper).

```bash
./run.sh programs/eta_vs_z.cc \
  -i data/jets/raw/jets_native.root \
  -o data/graphs/eta_vs_z.root
```

Overlay canvases (probably outdated):

```bash
./run.sh programs/overlay.cc \
  -i data/graphs/performance.root \
  -o data/graphs/performance_overlays.root

./run.sh programs/overlay.cc \
  -i data/graphs/controls.root \
  -o data/graphs/controls_overlays.root \
  --no-errors --same-algorithm
```

## Data contracts between stages

`cluster_jets` writes these trees containing jet info to an output:

- `LabFrameAntiktTruth`
- `LabFrameCentauroTruth`
- `BreitFrameAntiktTruth`
- `BreitFrameCentauroTruth`
- `LabFrameAntiktReco`
- `LabFrameCentauroReco`
- `BreitFrameAntiktReco`
- `BreitFrameCentauroReco`

Each tree entry stores one clustered jet with branches for event and jet content (see bottom of jet_tools/include/types.h). Jet kinematics are taken from clustered FastJet jets, and constituent IDs come from PseudoJet::user_index(), which is set to the original input-particle index and used to map back to the source particle collection.

`diff_frame/match_jets` writes:

- `best_matches_truth` - Closest (Boosted Centauro, lab anti-kt) pair within dR.
- `all_matches_truth`  - All lab anti-kt jets within dR of the boosted Centauro jet.
- `best_matches_reco`  - Reconstructed version.
- `all_matches_reco`   - Reconstructed version.

`same_frame/match_jets` writes:

- `antikt_matches` (`event`, `truth_index`, `reco_index`, `dR`)
- `centauro_matches` (`event`, `truth_index`, `reco_index`, `dR`)