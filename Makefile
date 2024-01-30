GF90 = gfortran

XFLAG = -fopenmp  -ffree-line-length-none -O3  -march=native  -fimplicit-none   -fwhole-file 

LIBS = -llapack

BIN = Thomson 

FFFLAGS = -Wall -g -msse4.2 -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant

FFLAGS = -Wall -Wno-unused -Wno-unused-dummy-argument -O2 

help:
	@echo  " -----------------------------------------\n to make the code you should do :  \n ----------------------------------------- \n make Thomson |  to compiled with gfortran \n -----------------------------------------";\
 

Thomson:
	$(GF90) src/*.f90 -o $(BIN) $(FFLAGS) $(LIBS)

    
clean:
	$(RM)  $(BIN)  animation.mp4 data.out figure.png test fort.4 fort.3 energy.dat data_frame.dat

Thomson_d: 
	$(GF90) src/*.f90 -o $(BIN) $(FFFLAGS) $(LIBS)
 
Thomson_mpi:
	mpif90 src/*.f90 -o $(BIN) $(FFLAGS) $(LIBS) 