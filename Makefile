all: quartic

CXXFLAGS=-g -Wall -std=c++17 -I/usr/include/eigen3

quartic: quartic.o curves.o
	g++ -o $@ $^ -lGL -lGLU -lglut

curves.o: curves.cc curves.hh vector.hh
