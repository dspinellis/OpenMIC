#
# Copyright 2011-2012 Diomidis Spinellis
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#

ifdef RELEASE
# Production flags
CXXFLAGS=-O3
else
# Debug build flags
CXXFLAGS=-g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_GLIBCXX_DEBUG -D_GLIBCXX_CONCEPT_CHECKS
# For testing
CXXFLAGS+=-std=gnu++0x -DTEST
endif

OBJ=mine.o ExtensiblePartition.o Point.o entropy.o Partition.o

all: mine.exe
	./mine Spellman-100.txt
	#./mine linear.txt

mine.exe: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@

mine.exe: mine.cpp ExtensiblePartition.h Partition.h Point.h debug.h entropy.h

ExtensiblePartition.o: ExtensiblePartition.cpp
Point.o: Point.cpp

entropy.o: entropy.cpp

gdb: mine.exe
	gdb mine.exe

MINE:
	java -jar MINE.jar linear.csv 0
