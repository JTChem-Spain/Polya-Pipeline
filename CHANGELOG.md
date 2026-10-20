# Changelog

All notable changes to Pólya Pipeline are recorded here.
The project follows [semantic versioning](https://semver.org/).

## [1.0.0] — 2026-10-20

First public version.

### Counting and enumeration
- Pólya counting of the derivatives obtained by replacing ring CH groups with one or more kinds of
  atoms, with an optional maximum number of substituents.
- Symmetry from the molecular graph, cross-checked against the symmetry of the 3D geometry.
- One representative per derivative, shown with its lowest IUPAC locants; patterns that give the
  same compound are merged and both counts are reported.
- IUPAC numbering of single rings, ring assemblies and common fused parents (naphthalene,
  quinoline, isoquinoline, azulene, indole, benzofuran, benzothiophene, anthracene, phenanthrene,
  pyrene); other ring systems are labelled as program numbering.
- Only aromatic CH are substituted when the molecule also has other CH (e.g. CH3); positions can
  be chosen with atom maps.

### Files
- `.xyz` geometries (ETKDG + MMFF, UFF fallback) and optional Gaussian 16 `.gjf` inputs with a
  configurable route line; CSV table of the results.

### Interface
- Graphical interface: live drawing with numbered positions, field validation with plain-language
  errors, pattern count before running, cancellable runs, sortable results table, drawing of the
  selected derivative, in-app user guide (F1) and keyboard shortcuts.
- Command-line interface with the same options (`python derivatives.py --help`).

### Validation
- `tests/control_test.py`: 32 molecules checked against Pólya's count, a brute-force count of
  distinct compounds, hand-calculated values, IUPAC locants and the output files.
