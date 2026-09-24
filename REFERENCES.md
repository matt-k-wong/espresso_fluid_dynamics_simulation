# References and attribution

Ideas and equations are cited here; citation is not permission to
redistribute third-party PDFs, figures, or datasets. All plots in this
repository are our own, generated from released example data or synthetic
fixtures. No author or company endorsement is claimed, and no product
branding is used beyond factual identification of the user's machine type
(early-production La Marzocco GS3 AV, reservoir-fed).

## Papers (inspiration / context, not contributors)

- Waszkiewicz, R., Myck, F., Białas, Ł., Puciata-Mroczynska, M.,
  Dzikowski, M., Szymczak, P., and Lisicki, M., "Under pressure:
  Poroelastic regulation of flow in espresso brewing," Phys. Fluids 38,
  063113 (2026). https://doi.org/10.1063/5.0319611 — preprint:
  https://arxiv.org/abs/2512.21528 (v2 consulted as accessible text;
  final journal version not independently inspected). Evolving-porosity
  hydraulic description and experimental evidence. Used as context for what
  our frozen-flow milestone lacks (no evolving porosity, no wetting/gas).
- Cameron, M. I., Morisco, D., Hofstetter, D., Uman, E., Wilkinson, J.,
  Kennedy, Z. C., Fontenot, S. A., Lee, W. T., Hendon, C. H., and
  Foster, J. M., "Systematically Improving Espresso: Insights from
  Mathematical Modeling and Experiment," Matter 2(3), 631–648 (2020).
  https://doi.org/10.1016/j.matt.2019.12.019 — flow inhomogeneity context;
  not proof of any universal pressure threshold.

Research-paper authors are inspirations/sources unless they actually
contributed code to this software.

## Marbling mathematics (separate workspace, credited here only)

- Jaffer, A., "Mathematical Marbling,"
  https://people.csail.mit.edu/jaffer/Marbling/ — mathematical generation
  of marbling designs, animations, and analysis (rake/serpentine/scallop
  patterns, Oseen flow, Lamb–Oseen vortex, pigment transport).
- Lu, S., Jaffer, A., Jin, X., Zhao, H., and Mao, X., "Mathematical
  Marbling," IEEE Computer Graphics and Applications 32(6), 26–35 (2012).
  https://doi.org/10.1109/MCG.2011.51
- The `pst-marble` CTAN package (Jaffer with J. Gilg and L. Manuel,
  2018–2019) and related arXiv notes (e.g. arXiv:1702.02106,
  arXiv:1810.04646) are Jaffer's work, cited as background only.
- The marbling workspace is separate and stays completely outside any
  espresso release candidate; no marbling code, figures, or data are
  bundled here. Jaffer bears no responsibility for this implementation.

## Machine guidance (factual, not a capability claim)

- La Marzocco guidance on pre-brew / pre-infusion and pressure manipulation
  (home.lamarzoccousa.com). Distinguishes AV pre-brew controls from MP paddle
  manipulation. Does not establish which functions/retrofits exist on this
  early unit; reservoir operation supplies no mains-line pressure. Verify the
  exact manual/settings before proposing pressure or pre-brew experiments.

## Software dependencies (actual licenses, venv 2026-09-23)

- Runtime: numpy 2.5.3 — BSD-3-Clause (bundled 0BSD/MIT/Zlib/CC0-1.0 parts
  per License-Expression); scipy 1.18.1 — BSD-3-Clause with bundled
  OpenBLAS (BSD-3-Clause), LAPACK (BSD-3-Clause-Open-MPI), GCC runtime
  library (GPL-3.0-or-later WITH GCC-exception-3.1), libquadmath
  (LGPL-2.1-or-later). SciPy bundles affect binary redistribution notices.
- Test-only: pytest 9.1.1 — MIT (plus transitive iniconfig MIT, packaging
  Apache-2.0, pluggy MIT, pygments BSD-2-Clause). Build: setuptools, wheel.
- No dependency PDFs, figures, or datasets are redistributed. Check the
  installed METADATA files for full bundled notices before any publication.

## License references (no license installed by this task)

- https://opensource.org/osd
- https://choosealicense.com/licenses/apache-2.0/
- https://www.copyright.gov/ai/
