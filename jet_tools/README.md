jet_calib helper modules

This folder is for code factored out of `scripts/jet_calib.cc`.
Keep `scripts/jet_calib.cc` as the entrypoint (main and wiring), and move
pure logic and I/O helpers here over time.

Contains calibration.cc, calibration.h, config.h, matching.cc, matching.h, root_io.cc, root_io.h, types.h, rho.cc, rho.h
