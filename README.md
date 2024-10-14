# Installation

Currently, there are only two supported building and development environments:

- [FNAL machines](#fnal-machines)
- [CNAF machine](#cnaf-machine)

## FNAL machines
FNAL machines are Almalinux9 and the supported package manager is [Spack](https://spack.io/). Documentation can be found [here](https://spack.readthedocs.io/en/latest/).

### Installation
To build and install `sandreco` using spack, download the installation [file](../../wiki/files/install-sandreco-with-spack.sh) and run it.

```console
curl -O https://raw.githubusercontent.com/wiki/DUNE/sandreco/files/install-sandreco-with-spack.sh
source install-sandreco-with-spack.sh <spack installation folder>
```

When the script ends to run, you are ready to use `sandreco`. 

### Development
To develop `sandreco`, first install `edepsim` downloading the installation [file](../../wiki/files/install-edepsim-with-spack.sh) and run it.

```console
curl -O https://raw.githubusercontent.com/wiki/DUNE/sandreco/files/install-edepsim-with-spack.sh
source install-edepsim-with-spack.sh <spack installation folder>
```

Then, build `sandreco` using the following commands:

```console
mkdir <installation path>
cd <installation path>
git clone https://github.com/DUNE/sandreco.git
cd build & mkdir build
spack load gcc@12.2.0
spack load root@6.28.12
spack load edepsim@3.2.0
cmake ../sandreco/ -DCMAKE_INSTALL_PREFIX=.. \
-DCMAKE_INSTALL_RPATH="$(spack find --paths edepsim@3.2.0 | grep "edepsim@3.2.0" | awk -F' ' '{print $2}')/lib:$(spack find --paths root@6.28.12 | grep "root@6.28.12" | awk -F' ' '{print $2}')/lib/root:${PWD}/../lib"
make -j8
make install
```

## CNAF machine
On the CNAF machine run, [CMake](https://cmake.org/) is used to build `sandreco` using the following commands:

```console
mkdir <installation path>
cd <installation path>
git clone https://github.com/DUNE/sandreco.git
sed -i "s:set(CMAKE_CXX_STANDARD 17):set(CMAKE_CXX_STANDARD 14):g" sandreco/CMakeLists.txt
cd build & mkdir build
source /opt/exp_software/neutrino/env.sh
cmake ../sandreco/ -DCMAKE_INSTALL_PREFIX=..
make -j8
make install
source ../setup.sh
```

## Installation with ups [DEPRECATED]

```console
$ VERSION="v01_00_00"
$ QUAL="e20:prof"
$ source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
$ mrb newDev -v $VERSION -q $QUAL
$ source localProducts_larsoft_*/setup
$ mrb g sandreco
$ cd $MRB_BUILDDIR
$ mrbsetenv
$ mrb b -j8
$ mrb i
```

for a specific tag or branch do `mrb g -t $TAG sandreco` or `mrb g -b $BRANCH sandreco`

## sandreco

The `sandreco` project provides six executables:
- **Digitize** will perform digitization,
- **SANDECALClustering** will clusterize the ECAL DAQ digit in clusters of reconstructed cells,  
- **Reconstruct** will reconstruct tracks in STT and clusters in ECAL
- **Analyze** will identify particles and assign them a momentum
- **FastCheck** will produce a lot of plots to check everything is ok
- **Display** displays events

The executables exploit two libraries:
- **libUtils.so** for utilities
- **libStruct.so** for i/o

# Run

### Digitize
- Create digits of STT and cells of calorimeter

```console
$ Digitize <MC file> <digit file>
```
### SANDECALClustering 
- Create clusters of ECAL reconstructed cells (output `<cluster file>` is not an input argument) 

```console
$ SANDECALClustering -d <digit file>
```

# Description

SANDECALClustering takes as input the digitized photodetector signals in `<digit file>` and produces clusters of reconstructed cells in the ECAL. The output `<cluster file>` has the following structure: 

# `TTree tCluster`
```
cluster
cluster.tid
cluster.x (y, z) 
cluster.t
cluster.e
cluster.ax (ay, az) #apex
cluster.sx (sy, sz) #direction
cluster.varx (vay, varz) #variance
cluster.reco_cells 
```
Each `reco_cell` object has the following structure:
```
int id;
double z;
double y;
double x;
double l;
int mod;
int lay;
double e;
double t; 
dg_ps ps1; #photodetector 1 digitized photo-signal
dg_ps ps2; #photodetector 2 digitized photo-signal

```

### Reconstruct
- Track find and fit of STT track
- Clustering of calorimeter cells

```console
$ Reconstruct <MC file> <digiti file> <reco file>
```

### Reconstruct using Drift Circles method
- reconstruct muon track fitting drift circles
- digitization included in the exhecutable
- reconstruct only muons in the fiducial volume
- reconstruct only tracks with at least 5+5 fired wires

```console
$ ./build/bin/ReconstructNLLmethod -edep <EDEP file> -wireinfo tests/wireinfo.txt -o <reco file> --hit_time --signal_propagation
```

### Analyze
- Evaluate parameters of particles
- Evaluate neutrino energy

```console
$ Analyze <MC file> <reco file>
```

### FastCheck
- Produce several plots to check everything is ok

```console
$ FastCheck <root file> <pdf file>
```

### Display
- Display an event

```console
$ Display <event number> <MC file> <input file> [show trajectories] [show fits] [show digits]
```

# Data format

The description of the data format can be found [here](../../wiki/Data-Model)


# Contribute

- The code format can be find [here](../../wiki/Code-Formatting)
- The developing scheme is described [HowToDevelop.pdf](https://baltig.infn.it/dune/sand-reco/-/wikis/uploads/8b897fb0ea753ef767b96312bdf9ccac/HowToDevelop.pdf)

# Support

For any communication, please refer to [DUNE-ND-SAND-SOFTWARE@fnal.gov](mailto:DUNE-ND-SAND-SOFTWARE@fnal.gov)

[![BUILD AND TEST SANDRECO](https://github.com/DUNE/sandreco/actions/workflows/build-and-test-sandreco.yml/badge.svg)](https://github.com/DUNE/sandreco/actions/workflows/build-and-test-sandreco.yml)
