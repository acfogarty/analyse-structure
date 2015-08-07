#CPPFLAGS=-g $(shell root-config --cflags)
#CPPFLAGS=-std=c++11 -g
CPPFLAGS=-std=c++11
#LDFLAGS=-g $(shell root-config --ldflags)
#LDLIBS=$(shell root-config --libs)
#LDFLAGS=-g

SRCS=main.cpp analyseStructure.cpp 
#OBJS=$(subst .cc,.o,$(SRCS))
OBJS=main.o analyseStructure.o

analyseStructure.exe: $(OBJS)
	g++ $(LDFLAGS) -o analyseStructure.exe main.o analyseStructure.o $(LDLIBS) 

analyseStructure.o: analyseStructure.cpp analyseStructure.h Vec3D.h
	g++ $(CPPFLAGS) -c analyseStructure.cpp

main.o: main.cpp analyseStructure.cpp analyseStructure.h
	g++ $(CPPFLAGS) -c main.cpp
