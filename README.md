<<<<<<< HEAD

AP library

-----------------
Install
-----------------

Compile 
	make
Modify 'AP.sh' (change to your path)
	APPATH=/home/xiaodongli/software/APLike
Add this to ~/.bashrc (change to your path)
	source /home/xiaodongli/software/APLike/AP.sh


-----------------
Usage
-----------------

Use 
	ifort main.f90 $APlm 
or equivalently
	ifort main.f90 -lAP -I/home/xiaodongli/software/APLike/mods
to compile with AP libarary

-----------------
Test
-----------------

	cd test
	ifort test.f90 $APlm
	./a.out

You will get (ignoring the information print at the beginning)

 	For Lambda CDM model with omegam =   0.2700000    
	  Chisq value of AP method =    70.8979493819427
	
-----------------
How to compile with cosmomc
-----------------

modify two places.

If you want to use APLike in, e.g., calclike.f90, add $$APlm in the options

$(OUTPUT_DIR)/calclike.o: calclike.f90 Makefile
        $(F90C) $(F90FLAGS) -c calclike.f90 $$APlm -o $(OUTPUT_DIR)/calclike.o


Also, add in the build of cosmomc EXE

cosmomc: directories camb $(OBJFILES)
        $(F90C) -o ../cosmomc $(OBJFILES) $(LINKFLAGS) $(F90FLAGS)  $$APlm
=======
# APLike
Likelihood for AP
>>>>>>> 69c87f997b86cfa65b525527796a6558528f8078
