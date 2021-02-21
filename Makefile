CC=clang-10
CXX=dpcpp
CFLAGS=-g -O3 -ffast-math -Wno-macro-redefined -Wall
CXXFLAGS=-std=c++20 -Wno-deprecated-anon-enum-enum-conversion
INCLUDE= -Isrc $(shell pkg-config --cflags lal lalsimulation gsl)
LIBS= $(shell pkg-config --libs lal lalsimulation gsl) -lm -fopenmp

openmp: parareal.o LALSimIMREOBNRv2.o src/write_waveform.c
	$(CC) $(CFLAGS) $(INCLUDE) parareal.o LALSimIMREOBNRv2.o src/write_waveform.c -o write_waveform $(LIBS)

oneapi: parareal_oneapi.o LALSimIMREOBNRv2.o src/write_waveform.c
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(INCLUDE) parareal_oneapi.o LALSimIMREOBNRv2.o src/write_waveform.c -o write_waveform $(LIBS)

LALSimIMREOBNRv2.o : src/LALSimIMREOBNRv2.c
	$(CC) $(CFLAGS) $(INCLUDE) -c src/LALSimIMREOBNRv2.c -fopenmp

parareal.o : src/parareal.h src/parareal.c
	$(CC) $(CFLAGS) $(INCLUDE) -c src/parareal.c -fopenmp

parareal_oneapi.o : src/parareal.h src/parareal_oneapi.cpp src/gamma.hpp src/derivative.hpp
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(INCLUDE) -c src/parareal_oneapi.cpp 

.PHONY : clean
clean :
	-rm parareal_oneapi.o parareal.o LALSimIMREOBNRv2.o write_waveform
