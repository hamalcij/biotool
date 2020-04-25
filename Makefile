# biotool
# Bioinformatics Toolbox
#
# Makefile
# Copyright (c) 2020 Hamalčík Jan
#
# Makefile for test/main.cpp
#

PROG = biotool
CXX = g++
CXXFLAGS = -std=c++17 -Wall
SRC = test/main.cpp library/quickhull/QuickHull.cpp $(wildcard library/*.cpp)
OBJ = $(SRC:.cpp=.o)
DEP = $(OBJ:.o=.d)

$(PROG) : $(OBJ)
	$(CXX) -o $@ $^

-include $(DEP)

%.d : %.cpp
	$(CPP) $(CXXFLAGS) $< -MM -MT $(@:.d=.o) >$@

.PHONY: clean
clean:
	rm -f $(OBJ) $(DEP) $(PROG)
