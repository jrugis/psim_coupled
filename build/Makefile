###########################################################################################
##
## target platforms: mahuika or iMac
##
## mahuika assumes:   CHECK VERSION???
##   ml gimkl/2017a
##   ml Eigen/3.3.4
##   ml SUNDIALS/4.0.0-gimkl-2017a
##
## notes 
##    source code formating: clang-format -i ./src/*
##    to execute on iMac: mpiexec -n 8 ./psim_coupled
##
###########################################################################################
## C++
CPP_DEFS = -std=c++11
CPP = mpicxx

###########################################################################################
## NOTE: uncomment block for ONE of debug, profile or optimised

## debug
#CPP_FLAGS =  $(CPP_DEFS) -m64 -Og -Wall -Wextra

## profile
#CPP_FLAGS =  $(CPP_DEFS) -m64 -O3 -pg 

## optimised
CPP_FLAGS =  $(CPP_DEFS) -m64 -O3 -march=native -mtune=native -Wall -Wpedantic -g
#CPP_FLAGS =  $(CPP_DEFS) -m64 -Og -Wextra -march=native -mtune=native -Wall -Wpedantic -g

###########################################################################################
## NOTE: uncomment lines for ONE platform
#CPP_INCLUDES = -I/usr/include  # mahuika 
CPP_INCLUDES = -I/usr/local/include/eigen3 # iMac

CPP_LIBS = -lsundials_cvodes

###########################################################################################

EXE = ../run/psim_coupled
#OBJS = main.o utils.o cAcinus.o cLumen.o cCVodes.o cCell_flow.o cCell_calcium.o cCellMesh.o
#OBJS = main.o utils.o cAcinus.o cLumen.o cLumenTree.o cCell_calcium.o cCellMesh.o cCell_flow.o
OBJS = main.o utils.o cAcinus.o cLumen.o cCell_calcium.o cCellMesh.o cCell_flow.o cCVode.o cLSODA.o LSODA.o
SRC = ./src

all: $(EXE)

main.o: $(SRC)/main.cpp $(SRC)/global_defs.hpp $(SRC)/cAcinus.hpp $(SRC)/cCell_calcium.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

utils.o: $(SRC)/utils.cpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cAcinus.o: $(SRC)/cAcinus.cpp $(SRC)/cAcinus.hpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp $(SRC)/cLumen.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cLumen.o: $(SRC)/cLumen.cpp $(SRC)/cLumen.hpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp 
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

#cLumenTree.o: $(SRC)/cLumenTree.cpp $(SRC)/cLumenTree.hpp
#	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

LSODA.o: $(SRC)/libsoda/LSODA.cpp $(SRC)/libsoda/LSODA.h $(SRC)/libsoda/helper.h
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cLSODA.o: $(SRC)/cLSODA.cpp $(SRC)/cLSODA.hpp $(SRC)/libsoda/LSODA.h $(SRC)/cCell_flow.hpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cCVode.o: $(SRC)/cCVode.cpp $(SRC)/cCVode.hpp $(SRC)/cCell_flow.hpp $(SRC)/global_defs.hpp $(SRC)/utils.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cCell_flow.o: $(SRC)/cCell_flow.cpp $(SRC)/cCell_flow.hpp $(SRC)/cCVode.hpp $(SRC)/cCell_calcium.hpp $(SRC)/cCellMesh.hpp $(SRC)/global_defs.hpp $(SRC)/cLSODA.hpp $(SRC)/utils.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cCell_calcium.o: $(SRC)/cCell_calcium.cpp $(SRC)/cCell_calcium.hpp $(SRC)/cCellMesh.hpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp $(SRC)/cCell_calcium.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cCellMesh.o: $(SRC)/cCellMesh.cpp $(SRC)/cCellMesh.hpp $(SRC)/cCell_calcium.hpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

$(EXE): $(OBJS)
	$(CPP) $(CPP_FLAGS) $(CPP_LIBS) -o $@ $+

clean:
	rm *.o
