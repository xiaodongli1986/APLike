
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
Compile
-----------------

In src/APfuns_002.f90 files, search for 
  	character(len=charlen), parameter :: covmatdir = '/home/xiaodongli/software/APLike/covmat_files/', &
	    chisqdir = '/home/xiaodongli/software/APLike/chisqs/'
	character(len=charlen), parameter :: syscorfiledir = '/home/xiaodongli/software/APLike/'
	character(len=charlen), parameter :: data2pcffiledir = '/home/xiaodongli/software/APLike/'
change the path to your APLike path

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

This shows computing AP likelihood for a set of LCDM parameters
	cd test
	ifort test.f90 $APlm -mkl
	./a.out


	

	
-----------------
Compile with cosmomc
-----------------

modify two places.

If you want to use APLike in, e.g., calclike.f90, add $$APlm in the options

$(OUTPUT_DIR)/calclike.o: calclike.f90 Makefile
        $(F90C) $(F90FLAGS) -c calclike.f90 $$APlm -o $(OUTPUT_DIR)/calclike.o


Also, add in the build of cosmomc EXE

cosmomc: directories camb $(OBJFILES)
        $(F90C) -o ../cosmomc $(OBJFILES) $(LINKFLAGS) $(F90FLAGS)  $$APlm


-----------------
Run using cosmomc
-----------------
... Zhenyu Zhang, tell us!...


-----------------
Play with different settings
-----------------
In AP_funs002.f90, you can find

	integer, parameter :: mubins(N1) = (/ 20,21,22,23,24,25 /)

Using larger mubin means better constraint, but also having the risk of larger systematics
Recommending using 15-20 for conservative, 20-25 for aggressive
