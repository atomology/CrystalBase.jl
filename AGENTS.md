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
- Formatting is Runic, run through pre-commit. Run it on the files you touched only,
  `pre-commit run --files <files>`, not `--all-files`.

## Testing
- Run the full package test suite from the repo root with `julia --project -e 'using Pkg; Pkg.test()'`.
- If you change docs content or public APIs, also verify the docs build with `julia --project=docs docs/make.jl`.
- Tests run through TestItemRunner: write new tests as `@testitem` blocks, not bare `@testset`.
- Doctests (`jldoctest`) run only in the full suite; a filtered `runtests.jl <file>` run skips them.
  To regenerate their output, temporarily uncomment `fix=true` in the `doctest` call of `test/runtests.jl`.
- In tests, call non-exported APIs with module qualification, for example `CrystalBase.coupled_axes(...)`.

## Conventions
- A lattice is a 3×3 matrix whose **columns** are the lattice vectors, in Å.
- A reciprocal lattice has the reciprocal lattice vectors as columns, in Å⁻¹,
  and includes the 2π factor: `recip_lattice = 2π * inv(lattice)'`.
- Fractional coordinates are relative to those columns: `cart = lattice * frac`.

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
- Type names capitalize each word as written in prose: `SpaceGroup`,
  `LittleGroup` (two words, as in Crystalline.jl, even though the functions
  are `spacegroup`, `littlegroup`), `Subpath` (one word), `KPath` (the
  symbol `k` plus "path", as in Brillouin.jl).
- `cart` / `frac` — Cartesian (Å, Å⁻¹) / fractional (lattice-relative)
  coordinates; every position-like name says which it is when both exist.
- Conversions between two representations of the same quantity are
  `a_to_b` (`cart_to_frac`); building a result from a parameter of a
  different kind is `x_from_y` (`kgrid_from_density`).
- Overload an existing function by argument type rather than coining an
  input-output noun pair: `atomic_symbol(14)` and `atomic_symbol("Fe1")`,
  not `label_symbol`.
- Internal names (unexported helpers, struct fields, function bodies) may be
  shorter: community-standard abbreviations are fine, and Unicode (`Γ`, `k₁`)
  is fine inside function bodies. Unicode never appears in exported names or
  file names.
- Atom identity uses three words, one per concept, never interchanged:
  - **element** — the chemical element, stored as its atomic number
    (`atomic_number`, `atom_numbers`) and written as its symbol
    (`atomic_symbol`, `atom_symbols`: `"Fe"`). Applied to an atom, "symbol"
    always means the element symbol; other symbols say what they are
    (`international_symbol`). "Element" itself stays in prose: identifiers
    name the representation (number or symbol), because the bare word also
    means a group element or a matrix element.
  - **label** — the per-site string: the element symbol, optionally refined
    by a suffix (`Fe1`, `Fe_up`, `X1`), as CIF's `_atom_site_label`, the
    `ATOMIC_POSITIONS` of Quantum ESPRESSO and the `atoms_frac` block of
    wannier90 name sites. Per-site data
    is always `atom_labels` (`labels` locally), never `species` or `names`;
    two atoms are compared by label ("the same label").
  - **species** — one distinct label, as a member of the set of a crystal's
    distinct labels (`unique_species`). It is the set-level noun, as in the
    literature and in QE's `ATOMIC_SPECIES`/`ntyp`/`ityp` and VASP's
    species lines: an index into that set, its size, and anything keyed by
    it (`Dict(species => ...)`, "per species") use it. Like QE and VASP,
    `Fe1` and `Fe2` are two species of one element. (pymatgen's per-site
    `site.species` is the one usage we do not follow.)

  `name` never refers to atoms.
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
- Run the checks under Testing before opening the PR.
- Keep changes focused; avoid unrelated refactors in the same PR.
- Summarize user-visible API changes in PR description and update README/docs examples when relevant.
- Confirm examples and snippets still run when changing user-facing API behavior.
