# Installation

Three installation and development options are supported:

- using [spack](#using-spack) (recommended on FNAL machines)
- using [cmake](#using-cmake) (recommended on CNAF machine)

## Using spack
[Spack](https://spack.io/) is a package manager for supercomputers, Linux, and macOS. Documentation can be found [here](https://spack.readthedocs.io/en/latest/).

### Installation
To install _sandreco_ using spack, download the installation [file](../../wiki/files/install-sand-with-spack.sh) and run it.

```console
curl -O https://raw.githubusercontent.com/wiki/DUNE/sandreco/files/install-sand-with-spack.sh
source install-sand-with-spack.sh <spack installation folder>
```

### Development
To develop _sandreco_ first install it as described [above](#installation). Then setup the spack development environment using the following commands:

```console
source <spack installation folder>/setup-env.sh
SANDENV="${SPACK_ROOT}/var/spack/environments/sand"
mkdir -p ${SANDENV}
spack env create -d ${SANDENV}
spacktivate ${SANDENV}
spack add sandreco
spack install
spack develop sandreco@01_00_00
spack concretize -f
spack install
```

The commands above should be run just once. The next time just use the following commands:

```console
source <spack installation folder>/setup-env.sh
SANDENV="${SPACK_ROOT}/var/spack/environments/sand"
spacktivate ${SANDENV}
```

Now you can develop the code. Once done, use the following commands to build _sandreco_:

```console
# change the code in ${SANDENV}/sandreco
spack install
```

## Using cmake
[CMake](https://cmake.org/) is the de-facto standard for building C++ code. It’s a powerful, comprehensive solution for managing the software build process. Documentation can be found [here](https://cmake.org/documentation/)

### Installation on FNAL machines
Assuming you have installed _edepsim_ in `<edepsim installation folder>`

```console
mkdir <installation path>
cd <installation path>
git clone https://github.com/DUNE/sandreco.git
mkdir build & cd build
spack load gcc@12.2.0
spack load root@6.28.12
cmake ../sandreco/ -DCMAKE_INSTALL_PREFIX=.. -DCMAKE_PREFIX_PATH=<edepsim installation folder>
make -j8
make install
```

#### Setup

```console
spack load gcc@12.2.0
spack load root@6.28.12
source <edepsim installation folder>/setup.sh
source <installation path>/setup.sh
```

### Installation on CNAF machine

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
```

#### Setup

```console
source /opt/exp_software/neutrino/env.sh
source <installation path>/setup.sh
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
