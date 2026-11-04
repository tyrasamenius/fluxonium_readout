# Change log

## 2026-11-04

- Removed the two-junction phase-basis methods and their dedicated operators
  from `Double_Junc_Fluxonium`; the notebook no longer includes cells that call
  the removed API.

## 2026-11-03

- Fixed `Ecs_from_caps_qir` in `Double_fluxonium_notebook.ipynb` to invert the
  full 3x3 capacitance matrix and derive diagonal charging energies and
  pairwise couplings from the inverse-matrix elements.
- Added the `_12_basis` suffix to the double-junction methods that operate in
  the two-junction phase basis, updated their notebook call sites, and added a
  matching device-resonator Hamiltonian builder.
- Moved the symmetric double-junction reduced-Hamiltonian helpers from
  `Fluxonium` to `Double_Junc_Fluxonium`. The latter now exposes the classical
  approximation, exact conditional BO correction, ground-band `A_00` and
  excited-band `A_a0` couplings, and the separate diagonal BH correction.
  `two_mode_h_symm` uses the exact conditional ground energy and leaves the
  BH correction out unless explicitly added.
- Added a reusable flux-sweep plotter that accepts a flux-dependent Hamiltonian
  factory, plus separate transition-sweep helpers for the BO and BH
  Hamiltonians.
- Extracted the conditional internal-mode Born-Oppenheimer analysis from
  `Double_fluxonium_notebook.ipynb` into a standalone Python module outside
  this repository. Its `plot_conditional_bo_analysis` function plots the
  lowest conditional bands, classical and geometrically corrected potentials,
  and neighboring-grid ground-state fidelity.

## 2026-10-28

- Updated `Double_Junc_Fluxonium.plot_wavefuncs_symm` in
  [`fluxonium_classes.py`](./fluxonium_classes.py) to plot the joint-basis
  probability `|psi(phi_q, n_i)|^2` instead of only the real part of each
  eigenvector. The color scale is shared across states, and the phase axis is
  labeled `phi_q`.
- Added `plot_wavefuncs_real_symm` to display signed real amplitudes with a
  diverging color scale centered at zero. Both plotting methods share the
  layout implementation and return their figure and axes.
- Added unit tests for the plotted data, labels, and color scales.
- Added this change log and workspace-level persistent instructions in
  [`.github/copilot-instructions.md`](./.github/copilot-instructions.md).

## Earlier conversation changes

- Added a phase-only `Fluxonium.get_two_mode_h` helper and routed the
  single-mode flux sweep through it. The helper has since been edited locally;
  consult the current source for its exact terms.
- Tried an absolute-value Fourier operator for `cos(Phi_i)` in the symmetric
  double-junction Hamiltonian, then removed that change at the user's request.
  The current operator is the signed charge-basis representation of
  `cos(Phi_i)`.

This log records assistant-authored changes identifiable from the conversation;
it does not attribute other pre-existing or local worktree edits.
