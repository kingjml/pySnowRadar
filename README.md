# pySnowRadar

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

 - [Batch-processing of multiple NSIDC L1b Deconvoluted SnowRadar products](./notebooks/batch_process_example.ipynb)
 - [Layer retrieval test of AWI SnowRadar product](./notebooks/retrieval_test_awi.ipynb)
 - [Layer retrievel test of OIB SnowRadar product](./notebooks/retrieval_test_oib.ipynb)

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
