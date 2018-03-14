#  ----------------------------
#  Make.inc for ATLAS3.10.2
#  ----------------------------

#  ----------------------------------
#  Make sure we get the correct shell
#  ----------------------------------
   SHELL = /bin/sh

#  -------------------------------------------------
#  Name indicating the platform to configure BLAS to
#  -------------------------------------------------
   ARCH = Corei264AVX

#  ----------------------------
#  Paths to various directories
#  ----------------------------
   BLDdir = /home/gerow/libatlas/ATLAS/midway_headnode
   SRCdir = /home/gerow/libatlas/ATLAS/midway_headnode/../
   INCAdir = $(BLDdir)/include
   INCSdir = $(SRCdir)/include
   BINdir = $(BLDdir)/bin
   LIBdir = $(BLDdir)/lib

   SYSdir = $(BLDdir)/tune/sysinfo
   GMMdir = $(BLDdir)/src/blas/gemm
   GMVdir = $(BLDdir)/src/blas/gemv
   GR1dir = $(BLDdir)/src/blas/ger
   L1Bdir = $(BLDdir)/src/blas/level1
   L2Bdir = $(BLDdir)/src/blas/level2
   L3Bdir = $(BLDdir)/src/blas/level3
   TSTdir = $(BLDdir)/src/testing
   AUXdir = $(BLDdir)/src/auxil
   CBLdir = $(BLDdir)/interfaces/blas/C/src
   FBLdir = $(BLDdir)/interfaces/blas/F77/src
   MMTdir = $(BLDdir)/tune/blas/gemm
   MVTdir = $(BLDdir)/tune/blas/gemv
   R1Tdir = $(BLDdir)/tune/blas/ger
   L1Tdir = $(BLDdir)/tune/blas/level1
   L3Tdir = $(BLDdir)/tune/blas/level3
   FLAdir = $(BLDdir)/src/lapack/reference
   ADdir  = $(SRCdir)/CONFIG/ARCHS

#  ---------------------------------------------------------------------
#  Name and location of scripts for running executables during tuning
#  ---------------------------------------------------------------------
   ATLRUN = $(BLDdir)/bin/ATLrun.sh
   ATLFWAIT = $(BLDdir)/bin/xatlas_waitfile

#  ---------------------
#  Libraries to be built
#  ---------------------
   ATLASlib = $(LIBdir)/libatlas.a
   CBLASlib = $(LIBdir)/libcblas.a
   F77BLASlib = $(LIBdir)/libf77blas.a
   LAPACKlib = $(LIBdir)/liblapack.a
   PTCBLASlib = $(LIBdir)/libptcblas.a
   PTF77BLASlib = $(LIBdir)/libptf77blas.a
   PTLAPACKlib = $(LIBdir)/libptlapack.a
   TESTlib = $(LIBdir)/libtstatlas.a

#  -------------------------------------------
#  Upper bound on largest cache size, in bytes
#  -------------------------------------------
   L2SIZE = -DL2SIZE=4194304

#  ---------------------------------------
#  Command setting up correct include path
#  ---------------------------------------
   INCLUDES = -I$(INCAdir) -I$(INCSdir) -I$(INCSdir)/contrib

#  -------------------------------------------
#  Defines for setting up F77/C interoperation
#  -------------------------------------------
   F2CDEFS = -DAdd_ -DF77_INTEGER=int -DStringSunStyle

#  ------------------------------
#  Architecture identifying flags
#  ------------------------------
   ARCHDEFS = -DATL_OS_Linux -DATL_ARCH_Corei2 -DATL_CPUMHZ=2600 -DATL_AVX -DATL_SSE3 -DATL_SSE2 -DATL_SSE1 -DATL_USE64BITS -DATL_GAS_x8664

TIDLIST=

#  -------------------------------------------------------------------
#  NM is the flag required to name a compiled object/executable
#  OJ is the flag required to compile to object rather than executable
#  These flags are used by all compilers.
#  -------------------------------------------------------------------
   NM = -o
   OJ = -c

#  ---------------------------------------------------------------------------
#  Various compilers and linkers.  The compilers recognized by ATLAS are:
#    F77 :  Compiler to compile f77 interface code with -- Should match what
#           user compile's his own code with; has no affect on ATLAS
#           performance
#    ICC :  Compiler to used to compile all non-kernels in ATLAS, including
#           C interface codes.  Should match whatever compiler user compiles
#           with.  Impact on performance should be minimal.
#    SMC :  Compiler to use to compile single precision matmul kernels
#    DMC :  Compiler to use to compile double precision matmul kernels
#    SKC :  Compiler to use to compile single precision non-matmul kernels
#    DKC :  Compiler to use to compile double precision non-matmul kernels
#    XCC :  Compiler to be used on the compile engine of a cross-compiler
#  These C compiler are mostly the same, but flags are often different.
#  Note that if you specify no F77 compiler, F77 will be set to a C compiler
#  to aid in linking.
#  ---------------------------------------------------------------------------
   NPROC=32
   CDEFS = $(L2SIZE) $(INCLUDES) $(F2CDEFS) $(ARCHDEFS) -m64 -DATL_DYLIBS -DATL_NCPU=$(NPROC)

   ICC = /usr/bin/gcc
   ICCFLAGS = $(CDEFS) -fomit-frame-pointer -mfpmath=sse -O2 -fno-schedule-insns2 -mavx -m64 -fPIC
   SMC = /usr/bin/gcc
   SMCFLAGS = -fomit-frame-pointer -mfpmath=sse -O2 -fno-schedule-insns2 -mavx -m64 -fPIC
   DMC = /usr/bin/gcc
   DMCFLAGS = -fomit-frame-pointer -mfpmath=sse -O2 -fno-schedule-insns2 -mavx -m64 -fPIC
   SKC = /usr/bin/gcc
   SKCFLAGS = -fomit-frame-pointer -mfpmath=sse -O2 -fno-schedule-insns2 -mavx -m64 -fPIC
   DKC = /usr/bin/gcc
   DKCFLAGS = -fomit-frame-pointer -mfpmath=sse -O2 -fno-schedule-insns2 -mavx -m64 -fPIC
   XCC = /usr/bin/gcc
   XCCFLAGS = $(CDEFS) -fomit-frame-pointer -mfpmath=sse -O2 -fno-schedule-insns2 -m64
   GCC = /usr/bin/gcc
   GCCFLAGS = -fomit-frame-pointer -mfpmath=sse -O2 -fno-schedule-insns2 -mavx -m64 -fPIC
   F77 = /usr/bin/gfortran
   F77FLAGS = -O -mavx -m64 -fPIC
   F77NOOPT = $(F77FLAGS) -O0   # turn off optimization
   SMAFLAGS =
   DMAFLAGS =
   CKC = $(SKC)
   ZKC = $(DKC)
   sKCFLAGS = $(CDEFS) $(SKCFLAGS)
   dKCFLAGS = $(CDEFS) $(DKCFLAGS)
   cKCFLAGS = $(CDEFS) $(SKCFLAGS)
   zKCFLAGS = $(CDEFS) $(DKCFLAGS)
   GOODGCC = /usr/bin/gcc -fomit-frame-pointer -mfpmath=sse -O2 -fno-schedule-insns2 -mavx -m64 -fPIC -m64 -fPIC
   KC = $(DKC)
   KCFLAGS = $(CDEFS) $(DKCFLAGS)
   LDFLAGS = -melf_x86_64
   F77SYSLIB = -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -lgfortran 
   BC = $(KC)
   NCFLAGS = $(KCFLAGS)

   CLINKER = $(KC)
   CLINKFLAGS = $(KCFLAGS)
   FLINKER = $(F77)
   FLINKFLAGS = $(F77FLAGS)
   FCLINKFLAGS = $(FLINKFLAGS)
   ARCHIVER = ar
   ARFLAGS  = r
   RANLIB   = echo

#  -------------------------------------
#  tar, gzip, gunzip, and parallel make
#  -------------------------------------
   TAR    = tar
   BZIP   = bzip2
   BUNZIP = bunzip2
   PMAKE  = $(MAKE) -j 32

#  ------------------------------------
#  Reference and system libraries
#  ------------------------------------
   FBLASlib = $(LIBdir)/libf77refblas.a
   FLAPACKlib = $(FLAdir)/lapack_$(ARCH).a
   SBLASlib = $(FBLASlib)  # should be serial sysblas
   BLASlib = $(FBLASlib)   # should be parallel sysblas
   SLAPACKlib =   # set to parallel system lapack
   SSLAPACKlib =  # set to serial system lapack
   LIBS = -lpthread -lm

#  --------------------------------------------------------------------
#  Compiler names for architectural defaults and flags to atlas_install
#  --------------------------------------------------------------------
   iccD = gcc
   smcD = gcc
   dmcD = gcc
   skcD = gcc
   dkcD = gcc
   xccD = gcc
   gccD = gcc
   f77D = gfortran
   INSTFLAGS = -1 0 -a 1 -l 1

#  -------------------------------------------------------------------
#  Dependence info for building optional external threading interfaces
#  -------------------------------------------------------------------
   sextthr =
   dextthr =
   cextthr =
   zextthr =
#  ---------------------------------------
#  Generic targets needed by all makefiles
#  ---------------------------------------
do_it: all
waitfile:
