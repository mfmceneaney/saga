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
cmake .. -DBUILD_DOXYGEN=FALSE
make
```
You should now have several executables in your `build` directory.

TODO: Add directions for adding dependency to another CMake project.

### Python3 Modules

For now just include these lines in your python code:
```python
import sys

sys.path.append('/path/to/saga/py')

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

## Getting Started
See the example configuration yamls for each C++ executable and the example python scripts in [tutorials/](tutorials/).

## Running an Analysis
Below is a detailed description of the project and how one might run an analysis.

### Design Philosophy
This library assumes your input data are stored in ROOT TTrees whose entries are event-level kinematics.

The provided C++ executables are entirely branch name agnostic.  These executables also assume configuration options will be passed in yaml file format which makes running analyses with different configuration option values extremely straightforward and easily organized.

Complementary to this aspect, the python libraries provide an easy interface for setting up job directories with different configurations and submitting jobs to slurm.  See, in particular, the `saga.orchestrate.create_jobs()` and `saga.orchestrate.submit_jobs()` methods.

The python libraries also provide some generic functionality for reading output asymmetry results from CSV files, aggregating results from job directories created by `saga.orchestrate.create_jobs()`, loading and computing systematics, and plotting asymmetry results and saving them to CSV.

### Examining Kinematics
The first step to any analysis should be to peruse your data.  You will  also need average kinematics in each bin if you plan to compare with theory results.  Use  the `getBinKinematics` executable (source code in [getBinKinematics.cpp](saga/src/getBinKinematics.cpp)) to get the statistics and average kinematics in each bin.

### Binning
Use the `findBinLimits` executable (source code in [findBinLimits.cpp](saga/src/findBinLimits.cpp)) to find the bin limits that will give roughly equal bin statistics for each binning variable in your dataset.

### Bin Migration
Detector reconstruction _smears_ reconstructed kinematics across the true generated kinematics.  To correct for this effect, one typically computes and then inverts a bin migration matrix.  The bin migration matrix is defined such that each entry

$M[i,j]=\frac{N[\text{Generated in bin }i\text{ and reconstructed in bin }j]}{N[\text{Generated in bin }i]}$.

Use the `getBinMigration` executable (source code in [getBinMigration.cpp](saga/src/getBinMigration.cpp)) to compute bin migration fractions and save the results to a CSV file.

### Invariant Mass Signal + Background Fits
Some physics channels isolate an asymmetry produced from a particle, such as a $\Lambda$ or $\pi^{0}$, which decays before hitting the detector array.  Since we only detect the final state particles one must examine the invariant mass spectrum of the 4-vector sum of the final state particles in order to isolate the contribution of the desired decay channel.  However, since a priori there is no way to know which particles we detect come from the same parent particle there will always be an irreducible combinatoric background under the signal distribution.

Only a $\Lambda$ invariant mass signal and background fit tuned specifically for CLAS12 RGA Fall 2018 Data is provided.  However, if you require a different signal and background fit the provided `analysis::applyLambdaMassFit()` method should be straightforward to expand upon.

#### Background Correction: Sideband Subtraction and sPlots

To correct for the background contributions to the final results there are two generally accepted methods.

The first, _Sideband Subtraction_ is very straightforward and simply uses a weighted subtraction of the asymmetry result computed in an arbitrary region(s) adjacent to the signal mass region, the sideband(s).  The weight $\varepsilon$ for the background is the relative fraction of background events within the signal region (also with arbitrary limits) computed from the invariant mass fit.  The assumption is that the asymmetry in the sideband region does not vary much from the asymmetry in the background under the signal peak, thus:

$A_{measured} = (1-\varepsilon)A_{signal} + \varepsilon A_{background}$.

Solving for $A_{signal}$:

$A_{signal} = \frac{A_{measured}-\varepsilon A_{background}}{(1-\varepsilon)}$.

The second method, $_sPlots$, which you may read about here: [arXiv:physics/0402083](https://arxiv.org/abs/physics/0402083), is a generalized form of sideband subtraction which computes event-level weights to produce the signal distribution in a variable which is _uncorrelated_ with the invariant mass variable.

### Fitting Asymmetries
Fitting an asymmetry $A$ may be done by minimizing a $\chi^{2}$ statistic or a Maximum Likelihood (ML) statistic.  In each case the acceptance effects of the detector must be corrected.

For a $\chi^{2}$ minimization fit (which is of necessity a binned fit) or a binned ML fit, this is done naturally by computing the binned asymmetry distribution

$A=\frac{N^{+}-N^{-}}{N^{+}+N^{-}}$

since the acceptance may reasonably be assumed to not depend on the helicity variable.

For an unbinned ML fit the acceptance naturally reduces to a relative luminosity factor between the positive and negative helicity subsets of the data ([H. Wollny, Thesis, University of Freiburg, 2010.](https://wwwcompass.cern.ch/compass/publications/theses/2010_phd_wollny.pdf), [G. Smith, Thesis, University of Glasgow, 2008.](https://theses.gla.ac.uk/5042/1/2013SmithPhD.pdf)).  This may usually be assumed to be $\simeq 1.0$.  Here the PDF takes the form:

$PDF(h,P_{b},\vec{x},\vec{a},\vec{d})=1+h\cdot P_{b} \cdot A(\vec{x},\vec{a},\vec{d})$

and $P_{b}$ is the polarization and $\vec{x}$, $\vec{a}$, $\vec{d}$ are the fit variables, asymmetry parameters, and depolarization variables (treated as independent variables).  In executables and functions provided by this project, the given asymmetry formula is converted internally to a PDF of this form and a simultaneous fit is done over the different helicity states.

Use the `getKinBinnedAsym` executable (source code in [getKinBinnedAsym.cpp](saga/src/getKinBinnedAsym.cpp)) to run a set of generically binned asymmetry fits and save the results to a CSV file.

### Injecting Asymmetries
To evaluate the effectiveness of your asymmetry extraction chain with all its fits and corrections, it is useful to _inject_ an artificial asymmetry into simulated data where the event-level truth is known for each event.  To do this one assigns a helicity $h \in [-1,1]$ for each event with a probability of a positive helicity taken from the asymmetry PDF.  In our provided executables this is done by generating a variable $r$ from a random uniform distribution for each event and then assigning a positive helicity if

$r<\frac{1}{2}PDF(h=1)=\frac{1}{2}(1+P_{b}\cdot Asymmetry(\vec{x}, \vec{a}, \vec{d}))$.

### Scaling MC for Future Experiments
To make projections of the uncertainties on planned asymmetry measurements for future experiments one must scale, in a sense "convert", uncertainties from simulated data to real data.  This requires one to compute the uncertainties on an existing simulation dataset and its corresponding real dataset.  These should provide a realistic comparison with the planned experiment.  For example, one might use CLAS12 RGC simulation and data to obtain CLAS12 RGH projections since both use a polarized target.

To scale these uncertainties, one first assumes Poissonian statistics so

$\delta N = \frac{1}{\sqrt{N}}$.

You must also compute ratio of the acceptance ratios

$R=\frac{N_{Reconstructed}}{N_{Generated}}$

for the two simulation samples and then the projected statistics in a given kinematic bin $i$ are given by:

$N_{New Data,i} = N_{Old Data,i} \cdot \frac{R_{New Sim.,i}}{R_{Old Sim.,i}} \cdot \frac{\mathcal{L_{New}}}{\mathcal{L_{Old}}}$,

where $\mathcal{L}$ denotes the integrated luminosity of each dataset.

In practice, rather than computing $N_{Generated}$ in each bin it is easier to divide $N_{Reconstructed}$ by the integrated cross-section since this should be directly proportional to $N_{Generated}$ and we take the ratio of the acceptance ratios so any extraneous units drop out anyway.

#

Contact: matthew.mceneaney@duke.edu
