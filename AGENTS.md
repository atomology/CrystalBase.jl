# AGENTS.md

## Scope
- These instructions apply to the entire repository.

## Development
- Use the latest stable Julia release.
- From the repo root, instantiate once with `julia --project -e 'using Pkg; Pkg.instantiate()'`.
- After dependency changes, keep the environment in sync with `julia --project -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'`.
- While iterating, run focused tests with `julia --project=test test/runtests.jl <test file name>...`.
  - Example: `julia --project=test test/runtests.jl type.jl`.
- This repo uses separate projects for `test` and `docs`; use the matching `--project` flag when needed.

## Testing
- Run the full package test suite from the repo root with `julia --project -e 'using Pkg; Pkg.test()'`.
- If you change docs content or public APIs, also verify the docs build with `julia --project=docs docs/make.jl`.
- In tests, call non-exported APIs with module qualification, for example `CrystalBase.coupled_axes(...)`.

## Naming
- Exported names use full words, no Unicode. Functions are verbs (with `!`
  for mutation) or the noun of what they return (`atomic_number`,
  `reciprocal_lattice`); types are nouns (`Crystal`, `KPath`).
- Separate words with `_` (`atom_labels`, `reciprocal_lattice`). Join words
  only when the pair is one established term, spelled as one word by the
  literature or by the Julia packages we interoperate with (Crystalline.jl,
  Spglib.jl, Base): `spacegroup`, `littlegroup`, `pointgroup`, `kpoints`,
  `kgrid`, `supercell`, `eigenvalues`. When unsure, use `_`.
- Prefixes and suffixes always take `_`, and a joined term stays whole:
  `n_atoms`, `n_kpoints`, `kpoints_cart`, `atom_positions_cart`.
  Predicates are the exception and join, as in Base (`isempty`, `haskey`).
- `cart` / `frac` — Cartesian (Å, Å⁻¹) / fractional (lattice-relative)
  coordinates; every position-like name says which it is when both exist.
- Conversions between two representations of the same quantity are
  `a_to_b` (`cart_to_frac`); building a result from a parameter of a
  different kind is `x_from_y` (`kgrid_from_density`).
- Overload an existing function by argument type rather than coining an
  input-output noun pair: `atomic_symbol(14)` and `atomic_symbol("Fe1")`,
  not `label_symbol`.
- Renames of exported names require prior approval.

## Code Style
- Use 4 spaces for indentation.
- Keep implementations small and type-stable where practical.
- Prefer explicit, readable linear algebra.
- Update docstrings and docs pages when changing public APIs in `src/`.
- Add or update tests for behavior changes, including edge cases for dimensions/units/conventions.
- If a local variable has the same name as a keyword argument, Julia lets you omit the keyword name in the call, for example `foo(x; y)`.

## PR checklist
- Recommended PR title format: `<short summary>`
- Ensure CI-equivalent checks pass locally:
  - Package tests
  - Docs build when docs/public APIs changed
- Keep changes focused; avoid unrelated refactors in the same PR.
- Summarize user-visible API changes in PR description and update README/docs examples when relevant.
- Confirm examples and snippets still run when changing user-facing API behavior.
