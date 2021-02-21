# GWSynth
Synthesizes audio from gravitational waveforms produced by binary black hole inspiral-merger-ringdown simulations.
The ODE solver stage of the simulation has been modified to use the Parareal parallel-in-time integration method: https://en.wikipedia.org/wiki/Parareal

## Dependencies
GWSynth uses the LAL and LALSimulation components of LALSuite: https://wiki.ligo.org/Computing/LALSuite

LAL and LALSimulation have the following dependencies:
* GSL
* fftw3
* hdf5

## Compilation
An OpenMP implementation of the code can be compiled with only a C compiler. 

There is also an Intel oneAPI DPC++ implementation that requires the oneAPI Base Toolkit: https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit.html 
