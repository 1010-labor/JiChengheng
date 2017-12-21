MPICXX=mpicxx
CFLAGS=-Wall -std=c++11

all:mpi_cfsfdp

mpi_cfsfdp:
	$(MPICXX) $(CFLAGS) -o mpi_cfsfdp mpi_cfsfdp.cpp

clean:
	rm -vf *.o mpi_cfsfdp
