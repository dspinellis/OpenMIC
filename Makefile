CXXFLAGS=-g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_GLIBCXX_DEBUG -D_GLIBCXX_CONCEPT_CHECKS
# For testing
CXXFLAGS+=-std=gnu++0x -DTEST

all: mine.exe
	./mine Spellman-300.txt xgrid ygrid

mine.exe: mine.cpp ExtensiblePartition.h Partition.h Point.h debug.h entropy.h
	$(CXX) $(CXXFLAGS) $< -o $@

gdb: mine.exe
	gdb mine.exe
