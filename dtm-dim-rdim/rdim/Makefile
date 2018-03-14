.SUFFIXES: .c .u .cpp
CC=g++ -w -Wunsafe-loop-optimizations

LIB=../lib
GSLWRAP_LIB=../gslwrap
LOCAL_DIR=../local
SPARSE_GSL_INC=../lib/include

CFLAGS = -DSPARSE -DMKL -O3 -DMKL_ILP64 \
	-DNDEBUG -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF  \
	-I ${SPARSE_GSL_INC} -I ${LIB} -I ${LIB}/math -I ${GSLWRAP_LIB}/include \
	-I ${GSLWRAP_LIB}/include/gslwrap -I ${LOCAL_DIR}/include \
	-I ${LIB}/util/gflags-1.1/src/gflags -I${LOCAL}/include \
        -I ${OPENBLAS_ROOT}/include \
	-I ${MKLROOT}/include \
	-m64 -fopenmp -frounding-math -msse4.2 -mfpmath=sse -march=native \
	-fno-trapping-math -fno-signaling-nans -fno-signed-zeros \
	-fgcse-after-reload -funsafe-loop-optimizations \
	-fgcse-sm -fgcse-las -fipa-pta -funswitch-loops 
        #-I ../opt-local/openblas/include \
	#-fprofile-use -fprofile-correction
	#-fprofile-correction -fprofile-generate -fprofile-arcs -fprofile-values

GDB_CFLAGS = -ggdb --enable-libunwind ${CFLAGS}

LDFLAGS = -Wl,-Map=ld_linkmap,-cref,-rpath -Wl,${LIB}/lib -L${LIB}/lib \
	-L${LOCAL_DIR}/lib -L${LOCAL_DIR}/lib/stl \
        -L../opt-local/openblas/lib -L/usr/lib/x86_64-linux-gnu/hdf5/serial \
	-lhdf5 -lhdf5_hl -lgsl -lopenblas -lm \
        -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_gnu_thread -lgomp -ldl -lpthread

RCC_LDFLAGS = -Wl,-Map=ld_linkmap,-cref,-rpath -Wl,${LIB}/lib -L${LIB}/lib \
	-L${LOCAL_DIR}/lib -L${LOCAL_DIR}/lib/stl \
	-L${OPENBLAS_ROOT}/lib \
	-lhdf5 -lhdf5_hl -lgsl -lopenblas -lm \
	-Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_gnu_thread -lgomp -ldl -lpthread
	#-L ../opt-rcc/openblas-rcc_amd/lib \

LOBJECTS = ss-lm.o gsl-wrappers.o data.o param.o util.o lda-seq.o lda.o params.o checkpoint.o main.o
POBJECTS = pardiso_solver.o
GFLAGS   = gflags.o gflags_reporting.o gflags_completions.o

all:	main solver

gflags:	${LIB}/util/gflags-1.1/src/gflags.cc
	$(CC) -c $(CFLAGS) \
	${LIB}/util/gflags-1.1/src/gflags.cc \
	${LIB}/util/gflags-1.1/src/gflags_reporting.cc \
	${LIB}/util/gflags-1.1/src/gflags_completions.cc

main:	$(LOBJECTS) gflags
	$(CC) $(CFLAGS) $(LOBJECTS) $(GFLAGS) -o main $(LDFLAGS)

rcc:	main_rcc solver_rcc

solver:	$(POBJECTS) gflags
	$(CC) $(CFLAGS)  $(POBJECTS) $(GFLAGS) -o pardiso_solver $(LDFLAGS)

solver_rcc:	$(POBJECTS) gflags
	$(CC) $(CFLAGS)  $(POBJECTS) $(GFLAGS) -o pardiso_solver $(RCC_LDFLAGS)

main_rcc:	$(LOBJECTS) gflags
	$(CC) $(CFLAGS)  $(LOBJECTS) $(GFLAGS) -o main $(RCC_LDFLAGS)

dbg:	$(LOBJECTS) gflags
	$(CC) $(GDB_CFLAGS) $(LOBJECTS) $(GFLAGS) -o main $(LDFLAGS)

clean:
	-rm -f *.o
	-rm -f main
