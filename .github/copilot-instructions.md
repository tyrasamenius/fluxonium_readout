# Persistent project preferences

Edit or add bullets below for guidance that should apply to future Copilot
requests in this repository.

- Treat Hamiltonian conventions, including signs and absolute values, as
  intentional modeling choices. Do not change them unless asked; explain the
  physical/model implications when relevant.
- Preserve unrelated local edits. Inspect the current implementation before
  changing it, and validate physics or numerical behavior with a focused check
  when practical.
- Keep reusable implementation and unit conversions in Python modules; keep
  experiment-specific parameters and analysis in notebooks.
- Keep variable names instructive and consistent with the physics context.
- When the user makes local edits, check that they align with the physics and give feedback without changing them. 
- Update the change log in `CHANGELOG.md` for any assistant-authored changes that are identifiable from the conversation.
