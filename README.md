# trimetric
Code & TeX archive for the paper "A variation on the Chamberlin trimetric map projection",
which introduces the new Matrix Trimetric map projection and compares it to the
Chamberlin Trimetric projection. See [here](https://github.com/brsr/trimetric/blob/main/tex/chamberlin_variation.pdf) for a PDF of the paper itself.

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

# Scripts

* `trimetric.py`: Generates figures needed for the TeX file. This seems to run faster in an iPython or Spyder terminal, probably because of matplotlib.
* `timing.py`: Times the two projections and compares to a baseline no-op.
* `invparams.py`: Used to check the convergence parameters of the inverse.
