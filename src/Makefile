
lib = ../lib/libAP.a

EXE1=AP_main
EXE2=AP_pp
EXE3=testfc
#EXE2 = AP_pp 
#OBJS = AP_funs.o
OBJS = AP_funs001.o


F90C = ifort -mkl
F90FLAGS = -heap-arrays

default: $(lib) $(EXE1) $(EXE2) $(EXE3)

$(lib): $(OBJS)
$(EXE1): $(OBJS) $(EXE1).o
$(EXE2): $(OBJS) $(EXE2).o

$(lib): $(OBJS)
	ar rc $(lib) $(OBJS)  
	cp *.mod ../mods/
$(EXE1): $(OBJS) $(EXE1).o
	$(F90C) -o $(EXE1) $(OBJS) $(EXE1).o 
$(EXE2): $(OBJS) $(EXE2).o
	$(F90C) -o $(EXE2) $(OBJS) $(EXE2).o 
$(EXE3): $(OBJS) $(EXE3).o
	$(F90C) -o $(EXE3) $(OBJS) $(EXE3).o 

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90 

#	$(F90C)  $(F90FLAGS)-c $*.f90 -lde -I/home/xiaodongli/software/delib/mods 

clean :
	rm -rf *.o *.mod *.a *.so  $(EXE1) $(EXE2) $(EXE3)
	rm -rf ../lib/*.a
	rm -rf ../mods/*.mod
	
