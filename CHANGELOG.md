# CHANGELOG

All notable changes to GENIE will be documented in this file.

## [Unreleased] - v1.2.2

### Documentation
- [DOC-001] (2025-12-10) Added `-s`/`--seed` option to README.md Parameters section. This option was functional but undocumented.
- [DOC-002] (2025-12-10) Added Memory Modes section to README.md explaining the three memory efficiency modes (0, 1, 2) and their configuration settings.
- [DOC-003] (2025-12-10) Added Config File Examples section with complete examples for Model G, G+GxE+NxE, and memory-efficient mode, plus a comprehensive parameters reference table.

### Bug Fixes
- [BUG-003] (2025-12-10) **Sample ID Matching**: GENIE now matches samples by FID_IID across all input files (genotype, phenotype, environment, covariate) and analyzes only the intersection. Previously, GENIE used position-based matching which could silently produce incorrect results if files had different sample orders. Added `--no-match-ids` flag and `no_match_ids` config option for legacy position-based matching (faster, requires pre-aligned files).

### Upgrades
- [UPG-001] (2025-12-10) Added validation for memory mode + model compatibility. GENIE now returns a clear error message instead of crashing (SIGSEGV) when memory-efficient mode is used with GxE models.
- [UPG-003] (2025-12-10) Added `-V`/`--version` flag and centralized version constant (`GENIE_VERSION`). Updated `--help` to document `-s`/`--seed` and `-V`/`--version` options.
- [UPG-004] (2025-12-10) Added unified `-mm`/`--memory-mode` flag and `memory_mode` config option. Simplifies memory mode configuration with values 0 (standard), 1 (memory-efficient), 2 (most memory-efficient).
- [UPG-007] (2025-12-10) Created centralized `validation.h` module for input validation. Provides clear warnings for unusual parameter values (num_vec, num_jack, nthreads) and validates configuration consistency.
