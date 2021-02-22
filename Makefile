CC=clang
CXX=dpcpp
CFLAGS=-g -fPIC -O3 -ffast-math -Wno-macro-redefined -Wall
CXXFLAGS=-std=c++20 -Wno-deprecated-anon-enum-enum-conversion
INCLUDE=-Isrc -Iinclude -Ilib $(shell pkg-config --cflags lal lalsimulation gsl)
LIBS=-llal -llalsimulation -lgsl -lm -fopenmp

openmp: libgwsynth.so src/write_waveform.c
	mkdir -p bin
	$(CC) $(CFLAGS) $(INCLUDE) src/write_waveform.c -o bin/write_waveform -Llib -lgwsynth -fopenmp

oneapi: libgwsynth_oneapi.so src/write_waveform.c
	mkdir -p bin
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(INCLUDE) src/write_waveform.c -o bin/write_waveform_oneapi -Llib -lgwsynth_oneapi -fopenmp -flto

libgwsynth_oneapi.so : parareal_oneapi.o LALSimIMREOBNRv2.o generate_waveform.o
	mkdir -p lib
	$(CXX) -fPIC -shared -Wl,-soname,libgwsynth_oneapi.so -o lib/libgwsynth_oneapi.so parareal_oneapi.o LALSimIMREOBNRv2.o generate_waveform.o $(LIBS)

libgwsynth.so : parareal.o LALSimIMREOBNRv2.o generate_waveform.o
	mkdir -p lib
	$(CC) -shared -Wl,-soname,libgwsynth.so -o lib/libgwsynth.so parareal.o LALSimIMREOBNRv2.o generate_waveform.o $(LIBS)

generate_waveform.o : src/generate_waveform.c
	$(CC) $(CFLAGS) $(INCLUDE) -c src/generate_waveform.c -fopenmp

LALSimIMREOBNRv2.o : src/LALSimIMREOBNRv2.c
	$(CC) $(CFLAGS) $(INCLUDE) -c src/LALSimIMREOBNRv2.c -fopenmp

parareal.o : src/parareal.h src/parareal.c
	$(CC) $(CFLAGS) $(INCLUDE) -c src/parareal.c -fopenmp

parareal_oneapi.o : src/parareal.h src/parareal_oneapi.cpp src/gamma.hpp src/derivative.hpp
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(INCLUDE) -c src/parareal_oneapi.cpp 

.PHONY : clean
clean :
	-rm -f generate_waveform.o LALSimIMREOBNRv2.o parareal.o lib/libgwsynth.so bin/write_waveform
