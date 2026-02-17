# centauro (new)

This directory contains the C++ analysis pipeline used to compare anti-k_t and Centauro jets on EDM4eic/EDM4hep event data. The workflow is:

1. Prepare events.
2. Cluster jets.
3. Match jets.
4. Produce analysis ROOT outputs and overlays.

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

From `new/`:

```bash
./build.sh
```

Useful options:

- `./build.sh clean`: remove `build/` before configure/build.
- `./build.sh fast`: skip configure, run build/install only.

Useful environment variables:

- `CENTAURO_BUILD_TYPE` (default `RelWithDebInfo`)
- `CENTAURO_WITH_EDM` (default `ON`)
- `CENTAURO_INSTALL_PREFIX` (default `build/install`)
- `BUILD_NPROC` (parallel jobs)

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

## End-to-end workflow

Run commands from `new/`.

### 1) Data prep

Reduce one raw EDM4eic file:

```bash
./run.sh programs/data_prep/reducer.cc \
  -i data/events/raw/events.root \
  -o data/events/filtered/events_skim.root
```

Batch download + reduce:

```bash
python3 programs/data_prep/get_batch_reco.py -l 0 -u 9 -c 5
```

Batch DIS/Breit cuts:

```bash
python3 programs/data_prep/cut_batch.py -l 0 -u 9
```

Merge reduced files:

```bash
./run.sh programs/data_prep/merge_data.cc \
  -o data/events/filtered/merged_pre.root \
  -i "data/events/filtered/*_skim.root"
```

Apply event-level DIS cuts:

```bash
./run.sh programs/data_prep/cut_events.cc \
  -i data/events/filtered/merged_pre.root \
  -o data/events/filtered/merged.root \
  --require-breit
```

### 2) Cluster jets

Native-threshold output (used by same-frame matching):

```bash
./run.sh programs/cluster_jets.cc \
  -i data/events/filtered/merged.root \
  -o data/jets/raw/jets_native.root \
  --threshold-frame native
```

Lab-threshold output (used by diff-frame matching):

```bash
./run.sh programs/cluster_jets.cc \
  -i data/events/filtered/merged.root \
  -o data/jets/raw/jets_lab.root \
  --threshold-frame lab
```

Defaults:

- `-R 0.6`
- `--max-eta 4`
- `--min-jet-pT 5`
- `--threshold-frame native`

### 3A) Same-frame analysis

Match truth/reco per algorithm:

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

### 3B) Diff-frame analysis

Cross-frame matching:

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

### 4) Additional plots

Eta vs z maps:

```bash
./run.sh programs/eta_vs_z.cc \
  -i data/jets/raw/jets_native.root \
  -o data/graphs/eta_vs_z.root \
  --min-jet-pT 2
```

Overlay canvases:

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

`cluster_jets` writes these trees to `data/jets/raw/*.root`:

- `LabFrameAntiktTruth`
- `LabFrameCentauroTruth`
- `BreitFrameAntiktTruth`
- `BreitFrameCentauroTruth`
- `LabFrameAntiktReco`
- `LabFrameCentauroReco`
- `BreitFrameAntiktReco`
- `BreitFrameCentauroReco`

`same_frame/match_jets` writes:

- `antikt_matches` (`event`, `truth_index`, `reco_index`, `dR`)
- `centauro_matches` (`event`, `truth_index`, `reco_index`, `dR`)

`diff_frame/match_jets` writes:

- `best_matches_truth`
- `all_matches_truth`
- `best_matches_reco`
- `all_matches_reco`

Main analysis outputs:

- `performance.root`: `summary`, `hists/`, `raw/`
- `controls.root`: `summary`, `hist/`, `antikt/`, `centauro/`
- `diff_frame_splitting.root`: `ratio/`, `eta_lab/`, `eta_breit/`, `nconst/`, `match/`, `shared_over_centauro/`
- `eta_vs_z.root`: `centauro/`, `antikt/`

## Pulling numeric summaries

Use:

```bash
python3 tools/read_outputs.py
```