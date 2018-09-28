# CYGNUS (work in progress)

This code is a collection of many years of work on directional detection culminating in a readout technology comparison for the CYGNUS project found in this paper [arXiv:1810.?????](https://arxiv.org/abs/1810.?????).

Contact me to complain about why the code won't work or whatever (ciaran.aj.ohare@gmail.com).

## Requirements

* python 2, and mpl_toolkits.basemap
* something to view the ipy notebooks
* if you wanna actually run the fortran code, then a fortran compiler. It should just work as is since it's not dependent on anything, just check the makefile for compiler information/flags etc.
* that's pretty much it, I've even given you an amputed minimiser code that I like

## Contents

* **code:** the main potatoes of the analysis, written in fortran (but hopefully in a friendly way), does most of the number-crunchy & MC stuff
* **data:** All the discovery reach curves get put here
* **neutrinos:** neutrino flux data
* **pixels:** lists of cartesian vectors from the HEALpix discretisation
* **plots:** all plots get put here
* **python:** plotting the outputs of the main code, as well as some copy-able functions in a more useful format
* **readouts:** readout technology performance data
