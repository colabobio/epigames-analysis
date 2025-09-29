# Analysis of epigames datasets

Analysis of the simulation data. Run the following Jupyter notebooks in order:

* 1-data-parsing
* 2-network-properties
* 3-super-spreader-analysis

The diffusion analysis needs be run in R:

* 4-risk-prediction

## Creating conda environment

The file requirements.txt list all the packages needed by these notebooks. It is recommended to use conda to create an environment with all this packages. 

First, install miniconda (or anaconda):

https://docs.anaconda.com/free/miniconda/

Clone this repo:

```
git clone https://github.com/colabobio/epigames-analysis.git
```

And the create the environment installing the listed requirements from the conda-forge channel:

```
cd epigames-analysis
conda create --name epigames --file requirements.txt --channel conda-forge
```
