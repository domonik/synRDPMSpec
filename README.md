# Minimal documentation

This repository exists for reproducibility of the corresponding publication.

*RAPDOR: Using Jensen-Shannon Distance for the computational analysis of complex proteomics datasets*



## System requirements

### Hardware requirements

This pipeline requires a standard computer with at least 16 GB of RAM. Howver 32 GB is reccomended

### Software requirements

#### OS
This pipeline requires a Linux system since some dependencies are not available for Windows or macOS it was tested using the following setup:

- Linux: Ubuntu 22.04 and Ubuntu 20.04

####  Dependencies

Dependencies are stored in the [environment.yml](./environment.yml) and in the yaml files in the [envs](./envs) directory. 
The specific versions used for the publication are stored in the `*.pin.txt` files


## Data preparation

You need the MaxQuant pre analyzed dataset in order to run the pipeline. Please contact the repository owner for required data.

First you need to download the following files and place them into a folder called Data

- https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.2c00759/suppl_file/pr2c00759_si_002.xlsx
- https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.2c00759/suppl_file/pr2c00759_si_003.xlsx

Note that you cannot download them programmatically consequently you need to use a browser to do so.


## How to run

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
