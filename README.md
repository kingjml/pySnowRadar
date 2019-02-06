

# pyWavelet (name TBC)

This is a package that Josh wants to be made legit.





## Usage and development instructions

### 1) Fork and/or Clone this repository

Do it.

### 2) Python environment creation

Get conda installed if you don't have it already:
```
# Windows x64 exe installer
https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe
# Linux x64 sh installer
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# MacOSX x64 pkg installer
https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.pkg
```

```
conda env create -f requirements.yml
conda activate py3-pyWavelet
```

### 3) Installing `pyWavelet`

```
cd /path/to/pywavelet/project/root/dir
conda activate py3-pyWavelet # skip this if you're coming from Step 2
pip install . 
```

Now you can import the `pyWavelet` module and its associated submodules, or checkout the jupyter notebooks within this repo.

### (Optional) Testing
#### Setup

```
cd /path/to/pywavelet/project/root/dir
```

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