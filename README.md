[![ROOT-latest](https://github.com/mfmceneaney/saga/actions/workflows/docker-image.yaml/badge.svg)](https://github.com/mfmceneaney/saga/actions/workflows/docker-image.yaml)

# Spin Asymmetry Generic Analysis (SAGA)

Run a generic SIDIS asymmetry analysis with CLAS12 data.

Use [CLAS12-Analysis](https://github.com/mfmceneaney/CLAS12-Analysis.git) or your own software to produce the input ROOT trees with event by event kinematics selecting all unique $e^{-}+X$ combinations.

## Install from source

Begin by cloning the repository:
```bash
git clone --recurse-submodules https://github.com/mfmceneaney/saga.git
```

### Prerequisites

* Python >=3.7.3
* A compiler with C++11 support
* Pip 10+ or CMake >= 3.4 (or 3.14+ on Windows, which was the first version to support VS 2019)
* [ROOT](https://root.cern.ch)

### C++ Libraries
This is a CMake project so you can build wherever but this is probably the simplest way to go.  From the top-level project directory run:
```bash
cmake -S . -B build -DBUILD_DOXYGEN=FALSE
cmake --build build
cmake --install build --prefix $PWD/bin
```

You should now have several executables in your `bin` directory.

You can manually source the environment scripts or add the following to your startup script to include these executables on your `$PATH`:
```bash
# Add SAGA environment variables: https://github.com/mfmceneaney/saga.git
cd /path/to/saga
source bin/env.sh
cd -
``` 

### Python3 Modules

You may install via pip using the path to the top-level repository directory:
```bash
pip install -e /path/to/saga
```
Then you can import the libraries in your python code with:
```python
import saga
```

## Install via Docker

You may also install and run the project as a [Docker](https://www.docker.com) container.

Begin as above by cloning the repository:
```bash
git clone --recurse-submodules https://github.com/mfmceneaney/saga.git
```

Then build the project with:
```bash
docker build -t saga-project /path/to/saga
```
and run it with:
```bash
docker run --rm -it saga-project
```
The `--rm` option tells docker to remove the container and its data once it is shut down.
To retain the container data though, you can mount a local directory (`src`) to a directory (`dst`)
inside the container with the following:
```bash
docker run --rm -it -v <src>:<dst> saga-project
```
Once you start the container you should have the following environment variables:
- `SAGA_HOME`
- `SAGA_BUILD`
- `SAGA_BIN`
Furthermore, the project executables in `$SAGA_BIN` should be available from your `$PATH`.

## Documentation
Check out the documentation page on [Read The Docs](https://saga.readthedocs.io/en/latest/)!

### Building the Documentation

Additional prerequisites for building the documentation:
* [doxygen](https://www.doxygen.nl)
* [sphinx](https://pypi.org/project/Sphinx/) (available with pip)
* [sphinx_rtd_theme](https://pypi.org/project/sphinx-rtd-theme/) (available with pip)
* [breathe](https://pypi.org/project/breathe/) (available with pip)
* [myst-parser](https://pypi.org/project/myst-parser/) (available with pip)

To build the documentation run cmake with the `BUILD_DOXYGEN` option set to `TRUE`:
```bash
cd build
cmake .. -DBUILD_DOXYGEN=TRUE
make
```

## Getting Started
See the example configuration yamls for each C++ executable and the example python scripts in the tutorials folder.

## Running an Analysis
Below is a detailed description of the project and how one might run an analysis.

### Design Philosophy
This library assumes your input data are stored in ROOT TTrees whose entries are event-level kinematics.

The provided C++ executables are entirely branch name agnostic.  These executables also assume configuration options will be passed in yaml file format which makes running analyses with different configuration option values extremely straightforward and easily organized.

Complementary to this aspect, the python libraries provide an easy interface for setting up job directories with different configurations and submitting jobs to slurm.  See, in particular, the `saga.orchestrate` module.

The python libraries also provide some generic functionality for reading output asymmetry results from CSV files, aggregating results from job directories created by `saga.orchestrate.create_jobs()`, loading and computing systematics, and plotting asymmetry results and saving them to CSV.  See, in particular, the `saga.aggregate` module.

### Examining Kinematics
The first step to any analysis should be to peruse your data.  Use the `getBinKinematicsTH1Ds` to produce 1D histograms of your kinematic variables in each bin.  You will  also need average kinematics in each bin if you plan to compare with theory results.  Use  the `getBinKinematics` executable to get the statistics and average kinematics in each bin.

### Binning
Use the `findBinLimits` executable to find the bin limits that will give roughly equal bin statistics for each binning variable in your dataset.

### Bin Migration
Detector reconstruction _smears_ reconstructed kinematics across the true generated kinematics.  To correct for this effect, one typically computes and then inverts a bin migration matrix.  The bin migration matrix is defined such that each entry

$M[i,j]=\frac{N[\text{Generated in bin }i\text{ and reconstructed in bin }j]}{N[\text{Generated in bin }i]}$.

Use the `getBinMigration` executable to compute bin migration fractions and save the results to a CSV file.

### Invariant Mass Signal + Background Fits
Some physics channels isolate an asymmetry produced from a particle, such as a $\Lambda$ or $\pi^{0}$, which decays before hitting the detector array.  Since we only detect the final state particles, one must examine the invariant mass spectrum of the 4-vector sum of the final state particles in order to isolate the contribution of the desired decay channel.  However, since _a priori_ there is no way to know which particles we detect come from the same parent particle, there will always be an irreducible combinatoric background under the signal distribution.

The method `saga::signal::fitMass()` allows one to apply a generic signal + background fit to an invariant mass spectrum in one or more variables.  Starting parameters and even PDF formulas for both signal and background may be loaded in a bin dependent way by specifying the `massfit_yamlfile_map` argument mapping bin scheme ids from `saga::analysis::getKinBinnedAsym()` to paths for YAML files containing the appropriate arguments to be loaded by `saga::signal::fitMass()`.

#### Background Correction

To correct for the background contributions to the final results there are two generally accepted methods.

##### Sideband Subtraction

The first, _Sideband Subtraction_ is very straightforward and simply uses a weighted subtraction of the asymmetry result computed in an arbitrary region(s) adjacent to the signal mass region, the sideband(s).  The weight $\varepsilon$ for the background is the relative fraction of background events within the signal region (also with arbitrary limits) computed from the invariant mass fit.  The assumption is that the asymmetry $A_{BG}$ in the sideband region does not vary much from the asymmetry in the background under the signal peak, thus the observed asymmetry is:

$A_{\text{Obs.}} = (1-\varepsilon)A_{SG} + \varepsilon A_{BG}$.

Solving for the signal asymmetry $A_{SG}$:

$A_{SG} = \frac{A_{\text{Obs.}}-\varepsilon A_{BG}}{(1-\varepsilon)}$.

Of course, the simplest way to apply this is to fit the asymmetry in the signal region and in the sideband region separately, and then find the true signal asymmetry $A_{SG}$ from the above equation.

However, oftentimes many asymmetry terms will need to be fit simultaneously, and it becomes more reliable to fit both the signal and sideband simultaneously.  This approach also allows one to treat $\varepsilon$ as function of the independent variables.  By fitting the mass spectrum in several different bins across the independent variables of the asymmetry, one may assign values for $\varepsilon$ that depend on the bin.  In this case the asymmetry to be fit is:

$A_{\text{Obs.}} = (1-\varepsilon(\vec{x}))A_{SG}(\vec{x}) + \varepsilon(\vec{x}) A_{BG}(\vec{x})$, for events in the signal region, and:

$A_{\text{Obs.}} = A_{BG}(\vec{x})$, for events in the sideband region.

One may then simply read off the fitted values for $A_{SG}$.

##### sPlots

The second method, $_sPlots$, which you may read about here: [arXiv:physics/0402083](https://arxiv.org/abs/physics/0402083), is a generalized form of sideband subtraction which computes event-level weights to produce the signal distribution in a variable which is _uncorrelated_ with the invariant mass variable.

### Fitting Asymmetries
Fitting an asymmetry $A$ may be done by minimizing a $\chi^{2}$ statistic or a Maximum Likelihood (ML) statistic.  In each case the acceptance effects of the detector must be corrected.

For a $\chi^{2}$ minimization fit (which is of necessity a binned fit) or a binned ML fit, this is done naturally by computing the binned asymmetry distribution

$A=\frac{N^{+}-N^{-}}{N^{+}+N^{-}}$

since the acceptance may reasonably be assumed to not depend on the beam helicity $\lambda$ or target spin $S$, whichever is the spin state of interest.

For an unbinned ML fit the acceptance naturally reduces to a relative luminosity factor between the positive and negative helicity subsets of the data ([H. Wollny, Thesis, University of Freiburg, 2010.](https://wwwcompass.cern.ch/compass/publications/theses/2010_phd_wollny.pdf), [G. Smith, Thesis, University of Glasgow, 2008.](https://theses.gla.ac.uk/5042/1/2013SmithPhD.pdf)).  This may usually be assumed to be $\simeq 1.0$.  Here the PDF takes the form:

$PDF(\lambda,S_{||},\phi_{S},\vec{x},\vec{a},\vec{d}) = 1 + A_{UU} + \lambda \cdot \overline{\lambda^2} \cdot A_{LU} + S_{||} \cdot \overline{S_{||}^2} \cdot D_{T} \cdot A_{UL}$
$+ \lambda \cdot S_{||} \cdot \overline{\lambda^2} \cdot \overline{S_{||}^2} \cdot D_{T} \cdot A_{LL} + \overline{S_{\perp}^2} \cdot D_{T} \cdot A_{UT}(\phi_{S})$
$+ \lambda \cdot \overline{\lambda^2} \cdot \overline{S_{\perp}^2} \cdot D_{T} \cdot A_{LT}(\phi_{S})$.

$\vec{x}$, $\vec{a}$, $\vec{d}$ are the fit variables, asymmetry parameters, and depolarization variables (treated as independent variables). $S_{||}$ is the longitudinal target spin vector parallel to the beam direction, and $S_{\perp}$ is the transverse target spin vector.  $\phi_{S}$ is the azimuthal angle of the target spin in the $\gamma^*N$ Center of Mass frame.  $\overline{\lambda^2}$ is the luminosity averaged beam polarization, and $\overline{S_{||}^2}$ and $\overline{S_{\perp}^2}$ are the luminosity averaged longitudinal and transverse target polarizations.  $D_{T}$ represents the target dilution factor.  In executables and functions provided by this project, the given asymmetry formula is converted internally to a PDF of this form and a simultaneous fit is done over the different spin states.  For a dataset of length $N$, the likelihood parameter used for parameter optimization is:

$\mathcal{L}(\vec{a}) = \prod_{i=1}^{N} PDF(\lambda_{i},S_{||,i},\phi_{S,i},\vec{x}_i,\vec{a}_i,\vec{d}_i)$.

An _extended_ ML Fit simply introduces the normalization factor $\mathcal{N}(\vec{a})$ as an optimization parameter assuming a Poissonian distribution so that the extended likelihood becomes:

$\mathcal{L}(\vec{a}) = \frac{\mathcal{N}(\vec{a})^Ne^{-\mathcal{N}(\vec{a})}}{N!}\prod_{i=1}^{N} PDF(\lambda_{i},S_{||,i},\phi_{S,i},\vec{x}_i,\vec{a}_i,\vec{d}_i)$.

Use the `getKinBinnedAsym` executable to run a set of generically binned ML asymmetry fits and save the results to a CSV file.

### Injecting Asymmetries
To evaluate the effectiveness of your asymmetry extraction chain with all its fits and corrections, it is useful to _inject_ an artificial asymmetry into simulated data where the event-level truth is known for each event.  The injection algorithm proceeds as follows.
For each event, a random number $r\in[0,1)$, beam helicity $\lambda\in(-1,0,1)$, and target spin $S\in(-1,0,1)$ are all randomly generated.
A non-zero $\lambda$ and $S$ are generated with probabilities taken from the beam and target polarizations respectively:

$P(\lambda \neq 0) = \overline{\lambda^2}$, and

$P(S \neq 0) = \overline{S^2}$.

The event weight $w$ is computed as:

$w = \frac{1}{2} (1 + A_{UU} + A_{UT}(\phi_{S}) + \lambda \cdot (A_{LU} + A_{LT}(\phi_{S})) + S_{||} A_{UL} + \lambda S_{||} A_{LL})$,

where $A_{UU}$, $A_{UT}$ are the unpolarized and (only) transverse target spin dependent asymmetry terms
and $A_{LU}$, $A_{LT}$, $A_{UL}$, and $A_{LL}$ are the asymmetry terms
dependent on beam helicity, target spin, or both.
Note that the asymmetry terms will taken from either the signal or background asymmetries
depending on whether the event has been marked as signal or background.
If $w>r$ the beam helicity and target spin values for that event are accepted,
otherwise all random values are regenerated and the process repeats until $w>r$.

### Scaling MC for Future Experiments
To make projections of the uncertainties on planned asymmetry measurements for future experiments one must scale, in a sense "convert", uncertainties from simulated data to real data.  This requires one to compute the uncertainties on an existing simulation dataset and its corresponding real dataset.  These should provide a realistic comparison with the planned experiment.  For example, one might use CLAS12 RGC simulation and data to obtain CLAS12 RGH projections since both use a polarized target.

To scale these uncertainties, one first assumes Poissonian statistics so

$\delta N = \frac{1}{\sqrt{N}}$.

You must also compute ratio of the acceptance ratios, that is, the ratio of reconstructed to generated counts:

$R=\frac{N_{Rec}}{N_{Gen}}$

for the two simulation samples and then the projected statistics in a given kinematic bin $i$ are given by:

$N_{New, Data,i} = N_{Old, Data,i} \cdot \frac{R_{New, Sim.,i}}{R_{Old, Sim.,i}} \cdot \frac{\mathcal{L_{New}}}{\mathcal{L_{Old}}}$,

where $\mathcal{L}$ denotes the integrated luminosity of each dataset.

In practice, rather than computing $N_{Gen}$ in each bin it is easier to divide $N_{Rec}$ by the integrated cross-section ($XS$) since this should be directly proportional to $N_{Gen}$ and we take the ratio of the acceptance ratios so any extraneous units drop out anyway.  In this case you would compute:

$N_{New, Data,i} = N_{Old, Data,i} \cdot \frac{N_{Rec, New, Sim.,i}}{N_{Rec, Old, Sim.,i}} \cdot \bigg{(}\frac{XS_{New, Sim.,i}}{XS_{Old, Sim.,i}}\bigg{)}^{-1} \cdot \frac{\mathcal{L_{New}}}{\mathcal{L_{Old}}}$.

#

Contact: matthew.mceneaney@duke.edu
