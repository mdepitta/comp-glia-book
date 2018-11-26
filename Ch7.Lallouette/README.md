# ChI network model simulation

The `src` directory contains the C++ code for simulating the network ChI model. The compiled code is used by the `RunSimulations.py` Python script. Figure-specific scripts like `Fig3.py` need to be called to plot the results. The following figures can be generated: 3C, 3D, 5C, 5D, 6E and 6F.

## Requirements

This software has been tested on Linux systems only.

The C++ code requires the `g++` compiler and the following packages:
- xutils-dev
- libgsl0-dev
- libboost-filesystem-dev
- libalglib-dev

The Python scripts all use Python 3 and require:
- [numpy](http://www.numpy.org/)
- [matplotlib](https://matplotlib.org/)

Make sure that these libraries are installed for Python 3 (not for Python 2.7).

## Compiling

```
cd src
make depend
make
```
Running `make depend` is usually only needed if you changed the C++ code.

## Running the simulations

```
python3 RunSimulations.py [-n NBCORES] [-f FIGNB...]
```
Without any arguments, it launches all simulations needed to generate all figures and uses all the available cores. The `-n` argument allows you to restrain the number of cores used while the `-f` argument is used to only run the simulations for a specific figure (3, 5, or 6).

The results are saved in the `data` folder; a subdirectory is created for each simulation groups (`data/Fig3Data` for Figure 3 for example) and simulations are distributed among 'jobs'. Each job then creates its own subdirectory (`data/Fig3Data/Sim_xxxxx`, with xxxxx being the random number generator seed).

`RunSimulations.py` contains some additional parameters like the number of repetitions of each network parameter combination or the number of simulations per job.

## Generating the figures

Once the corresponding simulations are done, the figures can be plotted using:
```
python3 Figure_3.py
python3 Figure_5.py
python3 Figure_6.py
```
These scripts take as an optional argument the path to the simulation group directory (e.g. `data/Fig3Data` for Figure 3).
