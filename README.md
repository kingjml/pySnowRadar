# pySnowRadar

![Build](https://github.com/kingjml/pySnowRadar/workflows/Build/badge.svg) ![Tests](https://github.com/kingjml/pySnowRadar/workflows/Tests/badge.svg) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4071801.svg)](https://doi.org/10.5281/zenodo.4071801)



A Python3 package to process data from CRESIS SnowRadar systems.

## From-Source Installation (PyPI and conda-forge packages coming soon)

Install and initialize conda:  https://docs.conda.io/en/latest/miniconda.html

Clone this repository, create the conda environment and install pySnowRadar:
  ```
  (base) $ git clone https://github.com/kingjml/pySnowRadar.git
  (base) $ cd ./pySnowRadar
  (base) $ conda env create -f exact_dev_env.yml
  (base) $ conda activate py3-pySnowRadar
  (py3-pySnowRadar) $ pip install . 
  ```

## Usage

Check out the Jupyter notebook examples for usage scenarios and code snippets:

 - [Batch-processing of multiple NSIDC L1b Deconvoluted SnowRadar products](https://github.com/kingjml/pySnowRadar/blob/master/notebooks/retrieval_batch_example.ipynb)
 - [Layer retrieval test of AWI SnowRadar product](https://github.com/kingjml/pySnowRadar/blob/master/notebooks/retrieval_awi_example.ipynb)
 - [Layer retrievel test of OIB SnowRadar product](https://github.com/kingjml/pySnowRadar/blob/master/notebooks/retrieval_oib_example.ipynb)
 
pySnowRadar does not validate interface or snow depth estimates. It is highly recommended that users compare outputs with measurements or references to quantify errors. Users should consider uncertainties including but not limited to surface roughness, salinity, and sidelobes.

For reference, the following papers (not an exhaustive list) describe uncertainties involved with the handling of SnowRadar data: [Panzer et al. (2013)]( https://www.cambridge.org/core/journals/journal-of-glaciology/article/an-ultrawideband-microwave-radar-for-measuring-snow-thickness-on-sea-ice-and-mapping-nearsurface-internal-layers-in-polar-firn/FF5CAF1FE9DBAD5847FF41995C1B1926), [Newman et al. (2014)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014JC010284), [Webster et al. (2014)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014JC009985), [King et al. (2015)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015GL066389), [Kwok et al. (2017)]( https://www.the-cryosphere.net/11/2571/2017/) 


## Development and Contributing

The following instructions are suitable for users who have already cloned this repository and would like to modify the inner workings of pySnowRadar.

From inside the local clone or fork, create a new branch where your modifications will reside:
  ```
  (py3-pySnowRadar) $ git checkout -b new_feature
  ```
  
After you make your modifications, you can test your changes by reinstalling pySnowRadar by using `pip` from within the local clone:

  ```
  # Make sure you're in the proper python environment!
  (py3-pySnowRadar) $ pip install . --upgrade
  ```
  
When you are satisfied with your changes, you may push your changes to github.com and open a pull request. For more information on pull requests, consult [Github's Documentation](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests)

## (Optional) Test-running
Test files are stored under the `tests` subdirectory and require the `pytest` and `coverage` packages.

#### Running tests

```
(py3-pySnowRadar) $ pytest 
```

#### Running tests and generating coverage reports
These commands run the tests and generate a coverage report for any untested files where the `Missing` columns indicates line numbers that still require testing

```
(py3-pySnowRadar) $ coverage run -m pytest
(py3-pySnowRadar) $ coverage report # to display the coverage report within a terminal
```

Sometimes it's nicer to see a html-ified coverage report, so use the following command:

```
(py3-pySnowRadar) $ coverage html 
```

This will generate a `htmlcov` folder containing the coverage report. Open `htmlcov/index.html` in a browser to see what code needs test coverage.
