# Spin Asymmetry Generic Analysis (SAGA)

Run a generic SIDIS asymmetry analysis with CLAS12 data.

Use [CLAS12-Analysis](https://github.com/mfmceneaney/CLAS12-Analysis.git) or your own software to produce the input ROOT trees with event by event kinematics selecting all unique $e^{-}+X$ combinations.

## Prerequisites

* Python >=3.7.3
* A compiler with C++11 support
* Pip 10+ or CMake >= 3.4 (or 3.14+ on Windows, which was the first version to support VS 2019)
* [ROOT](https://root.cern.ch)

## Installation

Begin by cloning this repository:
```bash
git clone --recurse-submodules https://github.com/mfmceneaney/saga.git
```

### C++ Libraries
This is a CMake project so you can build wherever but this is probably the simplest way to go:
```bash
mkdir build
cd build
cmake ..
make
```
You should now have several executables in your `build` directory.

TODO: Add directions for adding dependency to another CMake project.

### Python3 Modules

For now just include these lines in your python code:
```python
import sys

sys.path.append('/path/to/saga/py/saga')

from saga import orchestrate, aggregate
```

## Documentation
Check out the documentation page on [Read The Docs](https://saga.readthedocs.io/en/latest/)!

### Building the Documentation

Additional prerequisites for building the documentation:
* [doxygen](https://www.doxygen.nl)
* [sphinx](https://pypi.org/project/Sphinx/) (available with pip)
* [sphinx_rtd_theme](https://pypi.org/project/sphinx-rtd-theme/) (available with pip)
* [breathe](https://pypi.org/project/breathe/) (available with pip)

To build the documentation run cmake with the `BUILD_DOXYGEN` option set to `TRUE`:
```bash
cd build
cmake .. -DBUILD_DOXYGEN=TRUE
make
```

## Running an Analysis
TODO

### Design Philosophy

### Examining Kinematics
The first step to any analysis should be to peruse your data.

### Binning

### Invariant Mass Signal + Background Fits

### Background Correction: Sideband Subtraction and sPlots

### Fitting Asymmetries

### Injecting Asymmetries

### Orchestrating Jobs

### Aggregating Jobs

### Tabulating and Plotting Results

### Scaling MC for Future Experiments

#

Contact: matthew.mceneaney@duke.edu
