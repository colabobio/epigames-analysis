# OO WKU 2023 simulation

Analysis of the simulation data. Run the notebooks in order:

* 1-data-parsing
* 2-network-properties
* 3-behavioral-analysis
* 4-super-spreader-analysis
* 5-tensor-factorization
* 6-risk-prediction

Data can be downloaded from the following Zenodo repository:

https://zenodo.org/records/10674401

## Creating conda environment

The file requirements.txt list all the packages needed by these notebooks. It is recommended to use conda to create an environment with all this packages. 

First, install miniconda (or anaconda):

https://docs.anaconda.com/free/miniconda/

Clone this repo:

```
git clone https://github.com/colabobio/oo-wku23.git
```

And the create the environment installing the listed requirements from the conda-forge channel:

```
cd oo-wku23
conda create --name oo --file requirements.txt --channel conda-forge
```
