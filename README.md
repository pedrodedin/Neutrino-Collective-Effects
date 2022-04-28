# Neutrino Collective Effects Python Numerical Solver

[![DOI](https://zenodo.org/badge/363162761.svg)](https://zenodo.org/badge/latestdoi/363162761)


An Open Source code to workout numerical solutions for neutrino evolution considering neutrino-neutrino forward scattering (Collective Effects). Although it can be used for different scenarios, our focus is on the supernovae environment.

We implement the neutrino evolution equation in the formalism of polarization vectors. A better description of the neutrino evolution equation, numerical implementation, and first results, will be described in a future paper.

**Author**: Pedro Dedin Neto

For questions, comments, or bug reports, please contact Pedro Dedin Neto (pedroneto293@gmail.com).

## Download

You can download the soruce code or use the git clone command:

```bash
git clone https://github.com/pedrodedin/Neutrino-Collective-Effects.git
```

You also need to install the requested pyhton libraries ( [NumPy](http://www.numpy.org),  [SciPy](https://www.scipy.org), [Matplotlib](https://matplotlib.org/) ).

## Usage

For a quick guide use the Python notebook: [`Quick_Guide.ipynb`](Quick_Guide.ipynb)

In summary, the system of equations to be numerically solved are implemented in the files that begin with "ODE", each one corresponding to a different scenario (e.g. isotropic gas, bulb model emission). In each file, we have a function construct the initial conditions, another for the system of equations, and a final one to solve the system and give the final answer. Up to date, we have the following scenarios implemented:

* Isotropic and Monoenergetic Neutrino Gas:
* Isotropic Neutrino Gas with Spectral Distribution:
* **In Construction:** Bulb Model Emission

The file [`Auxiliar_Functions.py`](Auxiliar_Functions.py) contains auxiliary functions used in the numerical solver, such as matter potential profiles 	&lambda;(r), neutrino-neutrino profiles &mu;(r), functions to read the output of the numerical solver, etc.

The file [`Plots.py`](Plots.py) holds functions to do plots and the Python notebook [`Animations.ipynb`](Animations.ipynb) is used to create animations, such as the spectra evolution, and the polarization vectors precession.

## Notes
 
 This code is still a work in progress.
