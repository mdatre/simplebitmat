#CFLAGS= -O3 -std=gnu++0x -funroll-loops -m64 -lm -I.
CFLAGS= -O3 -std=c++11 -w -lm  -I.
#DEFINES= -DNUMTHREADS=`grep -c ^processor /proc/cpuinfo` #depends on underlying arch
#ASSERTS= -DNDEBUG
DEFINES= -DNUMTHREADS=8 #depends on underlying arch
#CFLAGSDBG=  -g3 -std=gnu++0x -m64 -lm -I.
CFLAGSDBG= -g -w -std=gnu++11 -m64 -lm -pthread -I.
BITMAT=./bin/bitmat
DICT=./bin/build_dictionary
DICTDBG=./bin/build_dictionary_dbg
TARGETDBG=./bin/bitmatdbg
CXX=nvcc

all: $(BITMAT)

$(BITMAT): main.o bitmatops.o 
	$(CXX) $(CFLAGS) -o $(BITMAT) bitmatops.o main.o 

main.o: main.cpp bitmat.hpp
	$(CXX) $(CFLAGS) -c main.cpp $(DEFINES)

bitmatops.o: bitmatops.cpp bitmat.hpp
	$(CXX) $(CFLAGS) -c bitmatops.cpp $(DEFINES)

bitmatdbg: main.cpp bitmatops.cpp bitmat.hpp 
	$(CXX) $(CFLAGSDBG) -o $(TARGETDBG) bitmatops.cpp main.cpp 

clean:
	rm $(BITMAT) $(TARGETDBG) $(DICT) *.o core
