# Optimal Transport
This repository contains the implementation of the algorithms described in my, Valentin Hartmann's, bachelor [thesis](https://arxiv.org/abs/1706.07403) and in the joint [paper](https://arxiv.org/abs/1706.07650) with Dominic Schuhmacher. That is the computation and visualization of the (semi-discrete) optimal transport for the Euclidean cost function from a continuous to a discrete measure.

In the folder [opt_transport](opt_transport) you can find the program that computes the optimal weight vector and the Wasserstein distance for two given measures whereas the folder [visualization](visualization) is dedicated to the additively weighted Voronoi diagrams. It comprises C++ code to generate the data necessary for drawing a Voronoi diagram for a set of points with associated weights and a MATLAB script which uses this data to carry out the job of the actual drawing.

The software was tested on Ubuntu 17.04 but it should work on most other Linux distributions without any problem.


## Installation
The programs make use of several software components that need to be installed first. The compilation and usage of the individual programs in this repository is explained in the READMEs in the subdirectories.

Open a shell window and install the following four packages by typing the corresponding command.

- CGAL for computing Voronoi diagrams
        sudo apt-get install libcgal-dev

- libLBFGS for minimizing the convex function Phi
        sudo apt-get install liblbfgs-dev

- CMake for compiling the sources
        sudo apt-get install cmake


## Typical workflow
To go from a source and a target measure to the Voronoi diagram inducing the optimal transport between them, we provide an example using the two sample measures provided in this repository. For details on how to use the individual programs see the corresponding READMEs.

First compile all programs as described. Then open a shell and enter the following:

        cd <root of this repository>
        ./opt_transport -l samples/mu.txt samples/nu.txt 3 ../weights.txt
        utilities/create_sites_file 3 ../samples/weights.txt ../samples/sites.txt
        visualization/create_diagram ../opt_transport/samples/sites.txt ../opt_transport/samples/intersections.txt

Change the paths in the [MATLAB script](visualization/plot_voronoi_diagram.m) to the sites and intersections file and run it.


## Thanks
- [CGAL](http://www.cgal.org/): A C++ library of compuational geometry algorithms.
- [libLBFGS](http://www.chokkan.org/software/liblbfgs/): A C implementation of the L-BFGS algorithm for minimizing convex functions.
