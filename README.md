This code-base is a fork of the dynamic topic models (DTM) and document influence models code (DIM).

This implementation fixes, optimizes and extends DTM / DIM provided by Gerrish & Blei (rDIM; see below). The code has been ported to use OpenMP an MPI-based dsitributed solver for rDIM (PARDISO). The code links against HDF5, Intel's (proprietary) MKL (w/ TPP) and OpenBLAS, in addition to GSL, linked to in the original. All of these should be **compiled from source** on the target system. On GNU systems, OpenBLAS must be compiled to use OpenMP (not pthreads) so that libgomp can correctly ration threads for BLAS calls.

Two changes to the code will affect DTM and DIM results compared to the older Gerrish & Blei DIM code. One is a bugfix where the original code initialized two RNGs, but only one of which respected the --rng_seed argument. Second, DTM and DIM used a Cholesky-based solver in their E-step. This has been changed to a more robust L-U Solver, which solves to higher precision. In most datasets exceeding ~1,000 documents this will produce values that differ at about 10e-8 precision.

This code also contains some optional optimizations of the original code, including fixing memory leaks, a checkpointing system, phi-update limiting in early iterations and the use of different BLAS libraries.

License is GPLv2. (see ./LICENSE), but see above for linking to required non-free libraries.

REQUIRES
MKL >=2017 (yes 2017) w/ TBB
GSL (current as of 2017-12-03)
OpenBLAS >=0.2.18
HDF5
64-bit CPU
A CPU with SSE 4.2+
libgomp or libiomp depending on compiler


INSTALLING:

When the stars are aligned

`make clean && make`

will create the `main` executable. `main --help` for help.

However ... it can be difficult to build and configure the model and its supporting libraries. The Makefile expects MKLROOT and OPENBLAS_ROOT to be set in the environment. At the time of writing, the required version of GSL (that supports sparse GMRES) was only available on GNU's Savannah git servers. Lastly, rDIM has a large linear solve in the model, which is accomplished by MKL's PARDISO MPI-based library. Other solvers work, though not well / fast / without tinkering with the make file.

This code has been tested primarily on Debian GNU/Linux and RHEL with icc/icpc (sorry) and to a less extent, gcc/g++. It has not been tested on MacOS, Unix, Windows or with clang, pgc or mpich.


Please cite as follows:

Dynamic Topics Models (DTM):

@inproceedings{blei2006dynamic,
  title={Dynamic topic models},
  author={Blei, David M and Lafferty, John D},
  booktitle={Proceedings of the 23rd international conference on Machine learning},
  pages={113--120},
  year={2006},
  organization={ACM}
}

Document Influence Model (DIM):

@inproceedings{gerrish2010language,
  title={A language-based approach to measuring scholarly impact},
  author={Gerrish, Sean and Blei, David M},
  booktitle={Proceedings of the 27th International Conference on Machine Learning (ICML-10)},
  pages={375--382},
  year={2010}
}

Regression-based Document Influence Model (rDIM):

@article{gerrish2010language,
  title={Measuring Discursive Influence Across Scholarship},
  author={Gerow, Aaron and Hu, Yuening and Boyd-Graber, Jordan and Blei, David and Evans, James},
  journal={Proceedings of the National Academy of Sciences},
  year={2018}
}

For questions about this code, contact Aaron Gerow (a.gerow@gold.ac.uk)
