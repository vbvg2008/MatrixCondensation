BIN = mc_p ge_p ge_scalapack
MODS=readMatrix.o


#module load OpenMPI
#module load ScaLAPACK

MODSRC=$(patsubst %.o,%.f,$(MODS))
MODMPI=$(patsubst %.o,%mpi.o,$(MODS))

DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
LIBS=-lscalapack -lblas

mc_p: mc_p.o $(MODMPI)
	mpifort -o $@ $^ $(LFLAGS)
ge_p: ge_p.o $(MODMPI)
	mpifort -o $@ $^ $(LFLAGS)
ge_scalapack: ge_scalapack.o $(MODMPI)
	mpifort -o $@ $^ $(LFLAGS)

mc_p.o: mc_p.f95 $(patsubst %.o,%mpi.o,$(MODS))
	mpifort -o $@ $< $(CFLAGS) -O3
ge_p.o: ge_p.f95 $(patsubst %.o,%mpi.o,$(MODS))
	mpifort -o $@ $< $(CFLAGS) -O3
ge_scalapack.o: ge_scalapack.f $(patsubst %.o,%mpi.o,$(MODS))
	mpifort -o $@ $< $(CFLAGS) -O3 $(LIBS)

$(MODS) : $(MODSRC)
	gfortran -o $@ $< $(CFLAGS)

$(MODMPI) : $(MODSRC)
	mpifort -o $@ $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f *.o *~ core $(BIN) *.mod
