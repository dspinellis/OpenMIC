CXXFLAGS=-g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_GLIBCXX_DEBUG -D_GLIBCXX_CONCEPT_CHECKS
# For testing
CXXFLAGS+=-std=gnu++0x -DTEST

OBJ=mine.o ExtensiblePartition.o Point.o entropy.o Partition.o

all: mine.exe
	./mine Spellman-300.txt xgrid ygrid

mine.exe: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@

mine.exe: mine.cpp ExtensiblePartition.h Partition.h Point.h debug.h entropy.h

ExtensiblePartition.o: ExtensiblePartition.cpp
Point.o: Point.cpp

entropy.o: entropy.cpp

gdb: mine.exe
	gdb mine.exe
