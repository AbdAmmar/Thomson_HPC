GF90 = gfortran

XFLAG = -fopenmp -llapack -ffree-line-length-none

BIN = Thomson

help:
	@echo  " -----------------------------------------\n to make the code you should do :  \n ----------------------------------------- \n make Thomson |  to compiled with gfortran \n ----------------------------------------- \n to make the input files for the code you should do : \n make Input   | \n -----------------------------------------"

Thomson:
	$(GF90) -c  include/*.f90 $(XFLAG)   $< ; \
  $(GF90) -c  $(BIN).f90 $(XFLAG)   $< ; \
  $(GF90) -o  Thomson *.o $(XFLAG)   $< ; \
  $(RM) *.o *.mod  $< 

Input:
	$(GF90) -o ING ING.f90  $<
    
clean:
	$(RM)  $(BIN) ING animation.mp4 data.out figure.png test fort.4 fort.3 energy.dat data_frame.dat
 