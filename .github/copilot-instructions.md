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

# Skills

### FIGURE
figure.figsize: 7, 5  

### LINES
lines.linewidth: 2    
lines.markersize: 4  

### LEGEND
legend.fontsize: 12
legend.frameon: False

### FONTS
font.size: 15

### AXES
axes.grid: False
axes.linewidth: 1
ytick.labelsize: 15
xtick.labelsize: 15



# Uses costum color cycle. Colors can be selected by color="C#", 
# where '#' is the number of the color in the cycle, e.g. 
# color="C0" for the first color in the cycle.

axes.prop_cycle: cycler('color', ['189BC4', 'FD8719', '762D34', '30A975', '021F52', 'D4002A', '70C6DD', 'FEBF7F', '9A5F65', '81CFAF', '5D7698', 'ED879D']) 

### GRID
grid.linewidth: 1
grid.linestyle: :
grid.color: k


### OTHER
savefig.dpi: 300
image.cmap: RdYlBu_r
