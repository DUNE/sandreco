# Description

SANDECALClustering takes as input the digitized photodetector signals in `<digit file>` and produces clusters of reconstructed cells in the ECAL. The output `<cluster file>` has the following structure: 

### `TTree tCluster`
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

# Installation

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


In the `bin` folder, there will be five executables:
- **Digitize** will perform digitization,
- **SANDECALClustering** will clusterize the ECAL DAQ digit in clusters of reconstructed cells, 
- **Reconstruct** will reconstruct tracks in STT and clusters in ECAL
- **Analyze** will identify particles and assign them a momentum
- **FastCheck** will produce a lot of plots to check everything is ok
- **Display** displays events

In the `lib` folder, there will be two libraries:
- **libUtils.so** for utilities
- **libStruct.so** for i/o

### Setup

```console
$ source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
$ setup edepsim v3_2_0b -q "e20:prof"
$ source setup.sh
```

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
