# pyWavelet

This is a package that Josh wants to be made legit.












## Environment creation

All installers are Python **3.7** at time of writing

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


## Testing
### Setup

```
cd /path/to/pywavelet/project/root/dir
```

### Running tests

```
pytest
```

### Generating test coverage report

```
coverage run --source pyWavelet -m pytest
coverage report
```