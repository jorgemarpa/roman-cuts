# Changelog

All notable changes to `roman-cuts` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-11-24

### Added
- Commissioning module (`roman_cuts.com`) for downloading and slicing Roman Commissioning (ROMCOM) data.
- MAST API query support (`query_mast_files`) with automated token resolution and product sorting.
- Multi-extension ASDF cube support for Roman data products.
- High-performance direct GWCS coordinate projections, fast bounding-box overlap slicing, and memory-mapped FITS/ASDF I/O.
- Parallel cutout extraction across multiple targets and SCAs with process and thread pooling.
- PEP 621 packaging with `uv` and `hatchling` backend.
- Unified dependency specifications across TRExS Roman toolchain.
- Standardized GitHub Actions CI workflows for pytest, ruff, mypy, and MkDocs deployment.

### Changed
- Migrated dependency management and packaging from Poetry to `uv`.
- Restructured documentation with MkDocs Material theme (`docs/README.md`, `docs/sim.md`, `docs/com.md`).
- Moved non-documentation notebooks and sandbox scripts out of the repository into `pipeline/`.

### Fixed
- Fixed WCS orientation and SCA coordinate order for simulated Roman mosaics.
- Corrected out-of-bounds padding and NaN handling for dithered cutout sequences.

## [0.3.1] - 2025-09-30

### Fixed
- Single FITS file cutout slicing bug and header metadata handling (#9).

## [0.3.0] - 2025-09-29

### Added
- Native ASDF cutout reader and writer support (#7).

## [0.2.2] - 2025-06-10

### Added
- Flux caching mechanism for accelerated re-slicing (#3).

## [0.2.0] - 2025-06-03

### Added
- Global WCS mapping across all 18 Roman WFI SCAs (#1).
- SCA coordinate transformation matrices and JSON coordinate lookup tables.

## [0.1.0] - 2025-03-31

### Added
- Initial release for creating image cutouts from `RImTimSim` simulations.
- Basic target coordinate extraction, multi-extension FITS support, and ASDF export.
