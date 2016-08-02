# **TopSpek** - A peak comparison and spectral decomposition tool)

Maintainer: Jonathan Williams


## Description

TopSpek is a peak comparison/decomposition tool for 1-D spectra, written in C.  It is mainly used for comparing simulated DSAM and RDM peak shapes with experimental data, though perhaps you can find other uses for it.  It compares peak shapes across simulated and experimental data in .mca (integer array) and .spe (radware) files.

The program works by minimizing the weighted chisq goodness of fit statistic between the input experiment and simulated dataset(s).  The minimization is performed by scaling the data and (if specified by the user) adding a background term (which may be constant, linear, or quadratic).  Multiple simulated datasets may be specified, and in this case the chisq of the sum of the datasets after scaling them individually will be minimized.

For a detailed example and description of how to specify parameters, see the included 'sample_parameters' file.


## Features

* Accepts .mca (integer array) and .spe (radware) files for experimental and simulated data.
* Simultaneously minimizes chisq goodness of fit statistic for one or multiple datasets with respect to the experimental data.
* Can add background of various forms and scale data to aid in chisq minimization.
* Scaling of data can be fixed to an absolute value or relative to other data, or be left as a free parameter.
* Plotting of fits is available via a built-in interface to gnuplot (program will still compile and run if gnuplot is not present, but plotting will be unavailable).
* Many other fitting and data processing options are availiable, see the included sample parameter file for details.


## Input Data Types

**.mca** - An .mca file is simply a 2D array of integers, with the first index denoting a spectrum number (up to 100) and the second index denoting a bin number (up to 32768).  See 'topspek_functions/save_data.c' for an example of how to write data to an .mca file.

**.spe** - An .spe file is the data type written by radware when using the 'ws' command in gf3.


## How to Install

Use 'make' to compile.  Optional data plotting requires gnuplot to be installed.

To run the program from anywhere, move the resulting 'topspek' executable to any directory under your $PATH environment variable.

Tested using gcc and GNU make on Ubuntu 14.04, Scientific Linux/CentOS 6, and Arch Linux (as of July 2016).  Should work on more or less any Linux distro.


## How to Run

The program runs from the command line, with only one argument: the path to a parameter file which sprecifies which files are to be fit and the parameters of the fit.

A sample parameter file is provided ('sample_parameters').


## Acknowledgments

This program uses the public domain gnuplot_i library by N. Devillard for displaying plots.
