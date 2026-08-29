# Versioning and File Roles

The source files originated as a chronological laboratory code series. Renaming them would break references in notes, datasets, and manuscripts, so the original `Photon_N[_variant].m` identifiers are preserved.

## Interpretation

- `Photon_N.m` is normally the primary simulation or study for case `N`.
- `Photon_N_1.m`, `Photon_N_2.m`, and later suffixes are follow-up analysis, aggregation, visualization, or parameter variants.
- `Photon_N_p.m` and `Photon_N_q.m` are named branches of the same case.
- `Photon_N_graphic.m` is a visualization-oriented variant.

Suffixes indicate lineage, not semantic software releases. Read the header of each file and the catalog before use.

## Preservation rule

Scientific logic and numerical parameters in historical files should not be changed without documenting the reason in Git history. Portable, reproducible replacements should be added as new curated workflows rather than overwriting the research record.

## Recommended future release convention

For publication-quality workflows, create a descriptive directory such as:

```text
workflows/localized_absorber_recovery/
  README.md
  run_simulation.m
  analyze_results.m
  config.m
  tests/
```

Tag validated repository releases using semantic versions such as `v1.0.0`.
