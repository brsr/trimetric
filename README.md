# trimetric
Code & TeX archive for the paper "A variation on the Chamberlin trimetric map projection",
which introduces the new Matrix Trimetric map projection and compares it to the
Chamberlin Trimetric projection. See [here](https://github.com/brsr/trimetric/blob/main/tex/chamberlin_variation.pdf) for a PDF of the paper itself.

Spun off from https://github.com/brsr/mapproj, which contains (among other
things) a Python prototype of the Matrix Trimetric Projection.

# Installation
This new projection is implemented in the framework of [PROJ](https://github.com/OSGeo/PROJ)

1. Download the source of PROJ version 8.
2. Copy `src/mattri.cpp` to `src/projections/` in PROJ.
3. Update these files to include `mattri` (just copy the `chamb` lines and change `chamb` to `mattri`):
* `src/Makefile.am`
* `src/pj_list.h`
* `src/lib_proj.cmake`
4. Compile and install PROJ as per its instructions.
5. Install [PyProj](https://github.com/pyproj4/pyproj) from source as per its instructions.

Check that the right version of PROJ is being recognized by PyProj with the command `pyproj.show_versions()`. You may need to remove other versions of PROJ and PyProj, or install into a fresh environment.

# Scripts
These scripts require Python 3.

* `trimetric.py`: Generates figures needed for the TeX file. This seems to run faster in an iPython or Spyder terminal, probably because of matplotlib.
* `timing.py`: Times the two projections and compares to a baseline no-op.
* `invparams.py`: Used to check the convergence parameters of the inverse.
