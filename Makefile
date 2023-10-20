all: quartic

CFLAGS=-g -Wall -std=c11 -pedantic
CXXFLAGS=-g -Wall -std=c++17 -pedantic -I/usr/include/eigen3 -DDEBUG -fsanitize=address

quartic: quartic.o curves.o nelder-mead.o reduction.o
	g++ -o $@ $^ -lasan -lGL -lGLU -lglut -lgmp

curves.o: curves.cc curves.hh vector.hh

nelder-mead.o: nelder-mead.cc nelder-mead.hh

reduction.o: reduction.c reduction.h
