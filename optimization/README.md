
In this assignment, we count distances between points in 3-dimensional space. One possibility to think of this is as cells, whose coordinates were obtained by diagnostics and which require further processing.


**Implementation**

Implement a program in C, called cell_distances, that uses OpenMP for parallelization and that:

Reads coordinates of cells from a text file “cells” in the current working directory.

Computes the distances between cells counting how often each distance occurs, approximating, i.e. rounding, truncating, or similar, to 2 decimal places.

Outputs to stdout a sorted list of distances with associated frequencies (including 0 fequencies at your will).
