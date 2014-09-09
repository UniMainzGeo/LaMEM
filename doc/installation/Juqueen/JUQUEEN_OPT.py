#!/usr/bin/env python
#
# this is taken from anton's earlier file, patched with some info found in the installed versions of PETSC on juqueen
# run with ./config/JUQUEEN_OPT.py
configure_options = [
  '--prefix=/homeb/hmz30/hmz300/LIB/petsc/3.4.2/JUQUEEN_OPT',
  'PETSC_ARCH=JUQUEEN_OPT',
  '--with-cc=mpixlc_r',
  '--with-cxx=mpixlcxx_r',
  '--with-fc=mpixlf77_r -qnosave',
  '--with-large-file-io',
# '--with-mpi-dir=/bgsys/drivers/ppcfloor/comm',  # required by BLACS to get mpif.h
# '--with-blas-lapack-lib=-L/soft/apps/LAPACK -llapack_bgp -L/soft/apps/LIBGOTO -lgoto',
  '--with-blas-lapack-lib=-L/bgsys/local/lib -lesslbg -L/bgsys/local/lapack/3.4.2/lib -llapack -L/bgsys/local/lib -lesslbg',
  '--download-superlu_dist=1',
  '--download-parmetis=1',
  '--download-metis=1',
# ML is deactivated here, as we will focuss on hypre instead, but see LaMEM doc on how to get ML working
#  '--download-ml=1',
#  '--download-ml=/homeb/hmz30/hmz300/ml-6.2.tar.gz',
  '--download-hypre=1',
  '--download-mumps=1',
  '--download-blacs=1',
  '--download-scalapack=1',
  '--with-x=0',
  '--with-is-color-value-type=short',
  '--with-shared-libraries=0',
  '-COPTFLAGS=-O3 -qstrict -qmaxmem=-1',
  '-CXXOPTFLAGS=-O3 -qstrict -qmaxmem=-1',
  '-FOPTFLAGS=-O3 -qstrict -qmaxmem=-1',
  '--with-debugging=0',
# autodetect on BGQ not working?
  '--with-fortran-kernels=1',
  '--with-batch=1',
  '--known-mpi-shared-libraries=0',
  '--known-memcmp-ok',
  '--known-sizeof-char=1',
  '--known-sizeof-void-p=8',
  '--known-sizeof-short=2',
  '--known-sizeof-int=4',
  '--known-sizeof-long=8',
  '--known-sizeof-size_t=8',
  '--known-sizeof-long-long=8',
  '--known-sizeof-float=4',
  '--known-sizeof-double=8',
  '--known-bits-per-byte=8',
  '--known-sizeof-MPI_Comm=4',
  '--known-sizeof-MPI_Fint=4',
#  '--known-mpi-long-double=1',
#  '--known-level1-dcache-assoc=0',
#  '--known-level1-dcache-linesize=32',
#  '--known-level1-dcache-size=32768',
#  '--download-hypre=1',
#  '--download-plapack=1',
#  '--download-parmetis=1',
#  '--download-umfpack=1',
#  '--download-triangle=1',
#  '--download-spooles=1',
#  '--download-superlu=1',
#  '--download-superlu_dist=1',
#  '--download-blacs=1',
#  '--download-scalapack=1',
#  '--download-mumps=1',
#  '--download-spai=1',
#  '--download-prometheus=1',
#  '--download-chaco=1',
  ]

if __name__ == '__main__':
  import sys,os
  sys.path.insert(0,os.path.abspath('config'))
  import configure
  configure.petsc_configure(configure_options)

# Extra options used for testing locally
test_options = []

