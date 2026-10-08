# IntraPhy 0.17 actual validation

- Source and installed wheel: 365 tests, 355 passed, 10 explicit missing-MAFFT skips.
- Runtime Python 3.13.5; NumPy 2.3.5, SciPy 1.17.0, Biopython 1.86 and NetworkX 3.6.1. The wheel used a separate environment with numeric dependencies exposed from the existing environment.
- Added 68 tests: 56 scope/candidate/visual contracts and 12 native-input extraction cases.
- CLI: parsimony, ER/ARD and foreground, 1 and 2 threads; 14 core result tables identical. 99 grouped SVGs parsed as XML.
- Four-taxon prepared-input gallery: 24 SVGs under all/high-coverage.
- Python 3.9 received syntax parsing only; other interpreters were not tested.
- No new real biological dataset run or statistical calibration. Missing external tools are not replaced with fake executables.

`validation.json` records scope, dependency versions and changed modules. Full logs accompany it. Older files in the parent validation directory are historical 0.16 evidence.
