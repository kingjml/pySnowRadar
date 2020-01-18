

<span style="color:red;font-size:  20pt">**!!! Alpha State - Do not use for publications !!!**</span>

# pySnowRadar

A Python3 package to process data from CRESIS SnowRadar systems. 

## Basic Installation

### For Windows Users

On Windows, there is usually some pain and suffering while managing the various geospatial dll linkages, so we recommend using [Anaconda](https://www.anaconda.com/distribution/) (specifically Miniconda) to install the required geospatial libraries *before* installing pySnowRadar:

Download and install Miniconda: 
- https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe

Either in the default `(base)` conda environment or a new environment, pre-install the geospatial library requirements:

 ```
 # this should pull in shapely, fiona and pyproj with proper linking
 conda install geopandas
 ```
 Then try installing pySnowRadar from PyPI which will bring with it any other package dependencies:

 ```
 pip install pySnowRadar
 ```

### For Linux Users

On Linux, the only requirement is you must have the Geospatial Data Abstraction Library (gdal) installed before using pySnowRadar. If you have root privileges, it is possible to install gdal at the system level using `apt`:

 ```
 sudo apt install gdal-bin
 ```

If you do not have root access, you may use `conda` to install gdal into the python environment:

 ```
 conda install gdal
 ```

After installing gdal, you may install pySnowRadar from PyPI which will automatically pull in the other package dependencies:

 ```
 pip install pySnowRadar
 ```

## Usage

Check out the Jupyter notebook examples for usage scenarios and code snippets:

 - [Batch-processing of multiple NSIDC L1b Deconvoluted SnowRadar products](https://github.com/kingjml/pySnowRadar/notebooks/batch_process_example.ipynb)
 - [Layer retrieval test of AWI SnowRadar product](https://github.com/kingjml/pySnowRadar/notebooks/retrieval_test_awi.ipynb)
 - [Layer retrievel test of OIB SnowRadar product](https://github.com/kingjml/pySnowRadar/blob/mike-dev/notebooks/retrieval_test_awi.ipynb)

## Development and Contributing

The following instructions are suitable for users who would like to modify the inner workings of pySnowRadar.

Clone this repository to your machine:
  ```
  git clone https://github.com/kingjml/pySnowRadar.git
  ```
Create a new branch where your modifications will reside:
  ```
  git checkout -b new_feature
  ```
After you make your modifications, you can reinstall pySnowRadar using your local clone by using `pip` from within the local clone of the github repository:

  ```
  # Make sure you're in the proper python environment!
  pip install . --upgrade
  ```

###  "Same-Env" installation

For convenience, we also provide `exact_dev_env.yml` that mirrors the `conda` environment used to develop pySnowRadar. If you encounter issues with the other installation steps, you may replicate our environment using the following commands using `conda` and `pip`:

```
(base) $ conda env create -f exact_dev_env.yml
(base) $ conda activate py3-pySnowRadar
(py3-pySnowRadar) $ pip install pySnowRadar
```

## (Optional) Test-running
Test files are stored under the `tests` subdirectory and require the `pytest` and `coverage` packages.

#### Running tests

```
pytest 
```

#### Running tests and generating coverage reports
These commands run the tests and generate a coverage report for any untested files where the `Missing` columns indicates line numbers that still require testing

```
coverage run -m pytest
coverage report # to display the coverage report within a terminal
```

Sometimes it's nicer to see a html-ified coverage report, so use the following command:

```
coverage html 
```

This will generate a `htmlcov` folder containing the coverage report. Open `htmlcov/index.html` in a browser to see what code needs test coverage.
