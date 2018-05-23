

EXE1 = test_om_chisqs 
EXE2 = test_om_w_chisqs
EXE3 = testfc_om_w_chisqs
EXE4 = testfc_CPL

F90C = ifort -mkl
F90FLAGS = -heap-arrays

default: $(EXE1) $(EXE2) $(EXE3) $(EXE4)


%: %.f90
	$(F90C) -o $*  $*.f90 $$APlm $(F90FLAGS) 


clean :
	rm -rf *.o *.mod *.a *.so  $(EXE1) $(EXE2) $(EXE3) $(EXE4)
	
