# scCRISPR-seq Perturbation Analysis Snakemake Workflow powered by Seurat's Mixscape
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow for performing perturbation analyses of pooled (multimodal) CRISPR screens with sc/snRNA-seq read-out (scCRISPR-seq) powered by the R package [Seurat's](https://satijalab.org/seurat/index.html) method [Mixscape](https://satijalab.org/seurat/articles/mixscape_vignette.html).

**If you use this workflow in a publication, don't forget to give credits to the authors by citing the URL of this (original) repository (and its DOI, see Zenodo badge above -> coming soon).**

![Workflow Rulegraph](./workflow/dags/rulegraph.svg)

Table of contents
----------------
  * [Authors](#authors)
  * [Software](#software)
  * [Methods](#methods)
  * [Features](#features)
  * [Usage](#usage)
  * [Configuration](#configuration)
  * [Examples](#examples)
  * [Links](#links)

# Authors
- [Stephan Reichl](https://github.com/sreichl)


# Software
This project wouldn't be possible without the following software and it's dependencies:

| Software       | Reference (DOI)                                   |
| :------------: | :-----------------------------------------------: |
| Mixscape       | https://doi.org/10.1038/s41588-021-00778-2        |
| Seurat         | https://doi.org/10.1016/j.cell.2021.04.048        |
| Snakemake      | https://doi.org/10.12688/f1000research.29032.2    |

# Methods
This is a template for the Methods section of a scientific publication and is intended to serve as a starting point. Only retain paragraphs relevant to your analysis. References [ref] to the respective publications are curated in the software table above. Versions (ver) have to be read out from the respective conda environment specifications (.yaml file) or post execution. Parameters that have to be adapted depending on the data or workflow configurations are denoted in squared brackets e.g. [X].

--- COMING SOON ---

# Features
The workflow performs all steps of the [Mixscape Vignette](https://satijalab.org/seurat/articles/mixscape_vignette.html) on all samples in the annotation file according to the parametrization in the config file.
- Calculation of local perturbation signatures
- Mixscape classification of perturbed cells versus cells with no detectable perturbation
- Visualization of Mixscape results
  - Perturbation scores of cells split by their mixscape classification
  - Posterior probability values in non perturbed and perturbed cells
  - if Antibody Capture was used: Violin plots of surface protein expression measurements split by perturbation classification of cells
- Analysis of perturbation responses with Linear Discriminant Analysis (LDA)
  -  Visualization using UMAP

# Usage
Read the [Mixscape Vignette](https://satijalab.org/seurat/articles/mixscape_vignette.html).

# Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

# Example
--- COMING SOON ---

# Links
- [GitHub Repository](https://github.com/epigen/mixscape_seurat/)
- [GitHub Page](https://epigen.github.io/mixscape_seurat/)
- [Zenodo Repository (coming soon)]()
- [Snakemake Workflow Catalog Entry](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/mixscape_seurat)
