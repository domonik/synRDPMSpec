# Data preparation

First you need to download the following files and place them into a folder called Data

- https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.2c00759/suppl_file/pr2c00759_si_002.xlsx
- https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.2c00759/suppl_file/pr2c00759_si_003.xlsx

Note that you cannot download them programmatically consequently you need to use a browser to do so.

# How to run

To run this pipeline you must install the required packages from the environment.yml file using conda:

```shell
conda env create -f environment.yml
conda activate synRDPMSpec
```

afterwards you can run the pipeline with the following command. Please adjust the cores parameter according to your system.


```shell
snakemake -s pipeline.smk --cores 4 --use-conda --configfile config.yml --rerun-incomplete
```

This will run the whole pipeline and place generated files in a Pipeline directory. 
