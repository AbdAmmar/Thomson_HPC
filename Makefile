FC = ifort #gfortran
LIBS = -mkl -qopenmp #-llapack -lblas

MPI_FLAGS = -DUSE_MPI

FFFLAGS = -Wall -g -avx512 -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant

#FFLAGS = -Wall -Wno-unused -Wno-unused-dummy-argument -O3 -fallow-argument-mismatch -g -mcpu=native -funroll-loops -cpp
FFLAGS = -O3 -march=native -fpp

SRC_DIR = src
INC_DIR = include
BLD_DIR = build
$(shell mkdir -p $(BLD_DIR))

F_SRC = $(wildcard $(SRC_DIR)/*.f90)
F_OBJ = $(F_SRC:$(SRC_DIR)/%.f90=$(BLD_DIR)/%.o)


BIN = Thomson

all: $(BIN)

$(BIN): $(F_OBJ)
	$(FC) $(FFLAGS) $(F_OBJ) -o $(BIN) $(LIBS)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(LIBS) -c -o $@ $<


Thomson_d: 
	$(FC) src/*.f90 -o $(BIN) $(FFFLAGS) $(LIBS)
 
Thomson_mpi:
	mpif90 src/*.f90 -o $(BIN) $(FFLAGS) $(MPI_FLAGS) $(LIBS) 

help:
	@echo  " -----------------------------------------"
	@echo  " to make the code you should do :         "
	@echo  " make Thomson                             "
 
.PHONY: clean
clean:
	$(RM) $(BIN) $(BLD_DIR)/* animation.mp4 data.out figure.png test fort.4 fort.3 energy.dat data_frame.dat Energy_parallel.out

