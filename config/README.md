# Configuration

You need one configuration file and one annotation file to run the complete workflow. If in doubt read the comments in the config and/or try the default values.

- project configuration (`config/config.yaml`): different for every project/dataset and configures the analyses to be performed.
- sample annotation (sample_annotation): CSV file consisting of two columns
    -  name: name of the dataset (tip: keep it short).
    -  data: absolute path to the Seurat object as .rds.

Set workflow-specific `resources` or command line arguments (CLI) in the workflow profile `workflow/profiles/default.config.yaml`, which supersedes global Snakemake profiles.