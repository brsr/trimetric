# trimetric
Code & TeX archive for the paper "A variation on the Chamberlin trimetric map projection",
which introduces the new Matrix Trimetric map projection and compares it to the
Chamberlin Trimetric projection.

Spun off from https://github.com/brsr/mapproj, which contains (among other
things) a Python prototype of the Matrix Trimetric Projection.

# Installation
This new projection is implemented in the framework of [PROJ](https://github.com/OSGeo/PROJ)
At some point I'll write an actual patch for this, but for now:

1. Download the source of PROJ, at least version 8. At the time of this writing, there's no release with that version, so just download the current version in development from the GitHub site.
2. Copy `src/mattri.cpp` to `src/projections/` in PROJ.
3. Update these files to include `mattri` (just copy the `chamb` lines).
* `src/Makefile.am`
* `src/pj_list.h`
* `src/lib_proj.cmake`
4. Compile and install PROJ as per its instructions.

# Commands

* `trimetric.py`: Generates figures needed for the TeX file. This seems to run faster in an iPython or Spyder terminal, probably because of matplotlib.
* `generate_pts_grid.py`: Generates a grid of points on the sphere.

For timing, where `degpts.txt` is the output of `generate_pts_grid.py`:
* `time proj +proj=mattri +lon_1=-80.0 +lon_2=-71.0 +lon_3=-35.0 +lat_1=9.0 +lat_2=-53.0 +lat_3=-6.0 +R=6371 < degpts.txt > /dev/null`
* `time proj +proj=chamb +lon_1=-80.0 +lon_3=-71.0 +lon_2=-35.0 +lat_1=9.0 +lat_3=-53.0 +lat_2=-6.0 +R=6371 < degpts.txt > /dev/null`
* `time cct +proj=noop < degpts.txt > /dev/null`
