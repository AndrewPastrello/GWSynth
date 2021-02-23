# GWSynth
Synthesizes audio from gravitational waveforms produced by binary black hole inspiral-merger-ringdown simulations.
The ODE solver stage of the simulation has been modified to use the Parareal parallel-in-time integration method: https://en.wikipedia.org/wiki/Parareal

## Dependencies
GWSynth uses the LAL and LALSimulation components of LALSuite: https://wiki.ligo.org/Computing/LALSuite

LAL and LALSimulation have the following dependencies:
* GSL
* fftw3
* hdf5 1.10

## Installation
An OpenMP implementation of the code can be compiled with only a C compiler as follows:

```
git clone https://github.com/AndrewPastrello/GWSynth.git --recurse-submodules
./build_lal.sh
```

Then make sure liblal and liblalsimulation are in your pkg-config and library paths and run:
```
make
```

This generates `bin/write_waveform`. For usage information run
```
write_waveform --help
```

There is also an Intel oneAPI DPC++ implementation that requires the oneAPI Base Toolkit: https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit.html

It can be compiled with:
```
make oneapi
```
