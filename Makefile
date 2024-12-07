GF90 = ifort #gfortran
LIBS = -mkl #-llapack -lblas

MPI_FLAGS = -DUSE_MPI

BIN = Thomson

FFFLAGS = -Wall -g -avx512 -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant

#FFLAGS = -Wall -Wno-unused -Wno-unused-dummy-argument -O3 -fallow-argument-mismatch -g -mcpu=native -funroll-loops -cpp
FFLAGS = -O3 -march=native -fpp

help:
	@echo  " -----------------------------------------\n to make the code you should do :  \n ----------------------------------------- \n make Thomson |  to compiled with gfortran \n -----------------------------------------";\
 

Thomson:
	$(GF90) src/*.f90 -o $(BIN) $(FFLAGS) $(LIBS)

    
clean:
	$(RM)  $(BIN)  animation.mp4 data.out figure.png test fort.4 fort.3 energy.dat data_frame.dat Energy_parallel.out

Thomson_d: 
	$(GF90) src/*.f90 -o $(BIN) $(FFFLAGS) $(LIBS)
 
Thomson_mpi:
	mpif90 src/*.f90 -o $(BIN) $(FFLAGS) $(MPI_FLAGS) $(LIBS) 
