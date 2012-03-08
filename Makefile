FLAGS = -L/usr/lib/X11
FLAGS += -L/usr/local/pgplot
FLAGS += -L$(MKLROOT)/lib/intel64
FLAGS += -I$(MKLROOT)/include/fftw
FLAGS += -O3
FLAGS += -fpe0
FLAGS += -march=corei7

FLAGS += -heap-arrays 4096
FLAGS += -no-prec-div
FLAGS += -xHost
FLAGS += -zero

#FLAGS += -parallel
#FLAGS += -par-num-threads=8

#FLAGS += -r8
#FLAGS += -axSSE4.1 
#FLAGS += -axSSSE3
#FLAGS += -fp-stack-check
#FLAGS += -fno-inline-functions
FLAGS += -traceback

LIBS = -lpgplot
LIBS += -lX11
LIBS += -lpng
LIBS += -lmkl_rt

COMP = ifort
all: Main Gauss

Main.o: Funcs.o TypesAndDefs.o Main.f90
	$(COMP) -c $(FLAGS) Main.f90 $(LIBS)
Gauss.o : Funcs.o TypesAndDefs.o Gauss.f90
	$(COMP) -c $(FLAGS) Gauss.f90 $(LIBS)
TypesAndDefs.o: TypesAndDefs.f90
	$(COMP) -c $(FLAGS) TypesAndDefs.f90 $(LIBS)
Funcs.o: TypesAndDefs.o Funcs.f90
	$(COMP) -c $(FLAGS) Funcs.f90 $(LIBS)
Main: Funcs.o TypesAndDefs.o Main.o
	$(COMP) $(FLAGS) -o Main Main.o Funcs.o TypesAndDefs.o $(LIBS)
Gauss: Funcs.o TypesAndDefs.o Gauss.o
	$(COMP) $(FLAGS) -o Gauss Gauss.o Funcs.o TypesAndDefs.o $(LIBS)
clean:
	rm *.mod *.o Main
remove:
	rm *.dat
