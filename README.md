# Installation

Currently, there are only two supported building and development environments:

- [FNAL machines](#fnal-machines)
- [CNAF machine](#cnaf-machine)

## FNAL machines
FNAL machines are Almalinux9 and the supported package manager is [Spack](https://spack.io/). Documentation can be found [here](https://spack.readthedocs.io/en/latest/).

### Installation
To install _sandreco_ using spack, download the installation [file](../../wiki/files/install-sandreco-with-spack.sh) and run it.

```console
curl -O https://raw.githubusercontent.com/wiki/DUNE/sandreco/files/install-sandreco-with-spack.sh
source install-sandreco-with-spack.sh <spack installation folder>
```

When script ends to run, you are ready to use _sandreco_. 

### Development
To develop _sandreco_, first install _edepsim_ downloading the installation [file](../../wiki/files/install-edepsim-with-spack.sh) and run it.

```console
curl -O https://raw.githubusercontent.com/wiki/DUNE/sandreco/files/install-edepsim-with-spack.sh
source install-edepsim-with-spack.sh <spack installation folder>
```

Then, install _sandreco_ using the following commands:

```console
mkdir <installation path>
cd <installation path>
git clone https://github.com/DUNE/sandreco.git
mkdir build & cd build
spack load gcc@12.2.0
spack load root@6.28.12
spack load edepsim@3.2.0
cmake ../sandreco/ -DCMAKE_INSTALL_PREFIX=..
make -j8
make install
source ../setup.sh
```

## CNAF machine
On the CNAF machine run, [CMake](https://cmake.org/) is used to build _sandreco_ using the following commands:

```console
mkdir <installation path>
cd <installation path>
git clone https://github.com/DUNE/sandreco.git
sed -i "s:set(CMAKE_CXX_STANDARD 17):set(CMAKE_CXX_STANDARD 14):g" sandreco/CMakeLists.txt
mkdir build & cd build
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

The `sandreco` project provides five executables:
- **Digitize** will perform digitization, 
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

### Reconstruct
- Track find and fit of STT track
- Clustering of calorimeter cells

```console
$ Reconstruct <MC file> <digiti file> <reco file>
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
