all: quartic

CXXFLAGS=-g -Wall -std=c++17 -I/usr/include/eigen3 -DDEBUG

quartic: quartic.o curves.o nelder-mead.o
	g++ -o $@ $^ -lGL -lGLU -lglut

curves.o: curves.cc curves.hh vector.hh

nelder-mead.o: nelder-mead.cc nelder-mead.hh
