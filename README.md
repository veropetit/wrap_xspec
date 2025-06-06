# wrap_xspec
Wrapper module to facilitate work with pyXspec

## Installation

Option 1:

You can clone the repository locally, and make sure that the wrap_xspec/wrap_xspec is in your PYTHON_PATH environment variable. 

The package is still in active development: please pull the repro often. 

Option 2:

You can use `pip install "git+https://github.com/veropetit/wrap_xspec`. 

The package is still in active development: please pull the repo often. 

## Requirements

This package assumes that you have xspec and pyXspec installed and acessible in your local environment. 

To verify:
* `xspec` in the command line should load xspec. 
* `python` and then `import xspec` in the command line should be working.
* `import xspec` in a jupyter notebook should be working. 

## Tutorials

Tutorials and their associated data are located in the Tutorial folder. If you are cloning the repository, I suggest that you make a copy of the Tutorial folder somewhere else to try them out (as to not clutter the repository and have to discard you changes when pulling the new version of the package). If you are using the pip method, simply download the Tutorial folder from the repo -- the folder contains all the necessary notebooks and example data. 