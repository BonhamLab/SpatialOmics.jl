# Tutorials

The tutorials are arranged around analysis decisions rather than package
implementation milestones. Start with the small, executable lessons and then
move to the technology-specific case studies.

## Core lessons

The following tutorials use synthetic data and run as part of the normal
Documenter build:

1. [Build a spatial dataset](@ref)
2. [Place an FOV in slide coordinates](@ref)
3. [Select acquisition sources and geometric ROIs](@ref)
4. [Assign transcripts and summarise expression](@ref)
5. [Persist a dataset safely](@ref)

These examples are deliberately small, but exercise the same public API used
for full experiments.

## Technology case studies

The Xenium and Visium tutorials use public datasets and committed rendered
figures. Their full inputs are several gigabytes, so routine documentation
builds display the curated outputs without downloading or recomputing them.
Each page also identifies the small native fixture used by tests and gives the
code used to regenerate it.

- [Xenium spatial transcriptomics](@ref)
- [Visium HD spatial transcriptomics](@ref)
- [Build a custom STARmap reader](@ref)

The [CosMx workflow](@ref) is a format guide rather than a collaborator-data
tutorial. Public-data figures can use the same pre-rendered approach when a
small redistributable fixture is selected.

This split keeps ordinary CI fast while keeping expensive examples
reproducible. Generated figures should record their public source dataset and
the code used to produce the committed asset; collaborator datasets are not
documentation inputs. Plotting dependencies for asset regeneration live in the
separate `docs/heavy` environment and are not installed by a normal docs build.
