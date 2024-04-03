# master-spreadsheet-datasets

A data pipeline used to generate a master spreadsheet of all the datasets shared by clusters in the Critical Zone Collaborative Network.

## How to use

This data analysis workflow uses Snakemake (installation instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)) as a pipelining tool to retrieve, clean and munge the spreadsheet to maximize readability.

First, create a Conda environment with all the required packages by running the following command: `
conda env create -f environment.yaml
`

Once in the new environment, we can execute the snakemake pipeline with this command: `snakemake --cores 1 -s Snakefile.smk --forceall`

When the jobs are done, the output master spreadsheet containing all cluster datasets will be in a newly created out folder in [3_munge/](3_munge/).

