# Requirements
- [ROOT](https://root.cern/)
- [edep-sim](https://github.com/ClarkMcGrew/edep-sim)

# Installation

### Get the code

```console
$ git clone git@github.com:DUNE/sandreco.git
```

### Build the binaries

```console
$ cd sandreco
$ mkdir build
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=./.. -DCMAKE_PREFIX_PATH=/path/to/edepsim/lib/cmake/EDepSim/ ./..
$ make
$ make install
```

In the `bin` folder, there will be five executables:
- **Digitize** will perform digitization, 
- **Reconstruct** will reconstruct tracks in STT and clusters in ECAL
- **Analyze** will identify particles and assign them a momentum
- **FastCheck** will produce a lot of plots to check everything is ok
- **Display** displays events

In the `lib` folder, there will be two libraries:
- **libUtils.so** for utilities
- **libStruct.so** for i/o

### Setup

```console
$ source setup.sh
```

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

<details><summary>Run production using batch system at CNAF</summary>

The batch system at CNAF is [HTCondor](https://htcondor.org/) (Documentation [here](https://htcondor.readthedocs.io/en/latest/)). To submit jobs, do:

```console
$ condor_submit -name sn-02.cr.cnaf.infn.it -spool submit.sub
```

The submit file is:

```
# Unix submit description file
# sleep.sub -- simple sleep job

executable              = submit.sh
arguments               = $(Item)
transfer_input_files    = /usr/lib64/libLHAPDF-6.2.1.so,macro.template.mac
log                     = $(ClusterId).$(Process).$(Item).log
output                  = $(ClusterId).$(Process).$(Item).out
error                   = $(ClusterId).$(Process).$(Item).err
should_transfer_files   = Yes
when_to_transfer_output = ON_EXIT
queue in 0, 21, 308
#queue from seq 300 399 | 
```

The instruction `queue in 0, 21, 308` will generate three jobs with $Item equals to 0, 21 and 308 respectively, while the instruction `queue from seq 300 400 |` will generate 100 jobs with $Item from 300 to 399. The variable $Item is then passed as [argument](https://baltig.infn.it/dune/nuev-generator/-/edit/master/README.md#L179) to the [script](https://baltig.infn.it/dune/nuev-generator/-/edit/master/README.md#L178) that will run on the node.

The script `submit.sh` is:

```
#!/bin/bash

# setup environment (i.e. ROOT and edep-sim)
source /opt/exp_software/neutrino/env.sh

# setup sand-reco
source /storage/gpfs_data/neutrino/users/mt/check_chain/sand-reco/setup.sh

# Output folder
# genie, edep-sim, digi and reco subfolders will
# be created and output files will be stored there
OUT=/storage/gpfs_data/neutrino/users/mt/check_chain/files

if [ ! -d ${OUT}/genie ]; then
  mkdir -p ${OUT}/genie
fi

if [ ! -d ${OUT}/edep-sim ]; then
  mkdir -p ${OUT}/edep-sim
fi

if [ ! -d ${OUT}/digi ]; then
  mkdir -p ${OUT}/digi
fi

if [ ! -d ${OUT}/reco ]; then
  mkdir -p ${OUT}/reco
fi

# NUEV-GENERATOR to be used
GEN=/storage/gpfs_data/neutrino/users/mt/check_chain/nuev-generator/bin/generate

# TOP VOLUME where neutrino interactions will be simulated
VOL=sand_inner_volume

# GDML geometry
GEO=/storage/gpfs_data/neutrino/SAND/GDML-GEO/SAND_opt2.gdml

# Neutrino FLUX
FLX=/storage/gpfs_data/neutrino/SAND/LBNF-NUBEAM-FLUX/
histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_LBNEND_fastmc.root

# GENIE tune
TUN=G18_02a_00_000

# XML cross-sections file
XSC=/storage/gpfs_data/neutrino/SAND/GENIE-XC/v3_00_06/NULL/G1802a00000-k250-e1000/data/gxspl-FNALbigger.xml

# Run number. It will be passed as argument
RUN=$1

# neutrino type
NUL=14

# number of simulated events (per job)
NEV=100

# Output file prefix
PREFIX=100-numu-InnerVol

# Output GENIE file prefix
# and edep-sim, Digitize, Reconstruct 
# output files
GENIE_OUT=${OUT}/genie/${PREFIX}
EDEP_OUT=${OUT}/edep-sim/${PREFIX}.${RUN}.edep-sim.root
DIGI_OUT=${OUT}/digi/${PREFIX}.${RUN}.digi.root
RECO_OUT=${OUT}/reco/${PREFIX}.${RUN}.reco.root

## NUEV-GENERATOR
${GEN} \
  -f histo:${FLX} \
  -g ${GEO} \
  -n ${NEV} \
  -o ${GENIE_OUT} \
  -p ${NUL} \
  -r ${RUN} \
  -t ${VOL} \
  --cross-sections ${XSC} \
  --tune ${TUN}

## GNTPC
if [ -f ${GENIE_OUT}.${RUN}.ghep.root ]; then
  gntpc \
    -f t2k_rootracker \
    -i ${GENIE_OUT}.${RUN}.ghep.root \
    -o ${GENIE_OUT}.${RUN}.gtrac.root
fi

## PREPARE MACRO FOR EDEP-SIM
sed -e "s:__NEV__:${NEV}:g" -e "s:__INPUT__:${GENIE_OUT}.${RUN}.gtrac.root:g" macro.template.mac > macro.${RUN}.mac

## EDEP-SIM
if [ -f ${GENIE_OUT}.${RUN}.gtrac.root ]; then
  edep-sim -C \
    -g ${GEO} \
    -o ${EDEP_OUT} \
    macro.${RUN}.mac
fi

## DIGITIZATION
if [ -f ${EDEP_OUT} ]; then
  Digitize \
    ${EDEP_OUT} \
    ${DIGI_OUT}
fi

## RECONSTRUCTION
if [ -f ${DIGI_OUT} ]; then
  Reconstruct \
    ${EDEP_OUT} \
    ${DIGI_OUT} \
    ${RECO_OUT}
fi

## ANALISYS
if [ -f ${RECO_OUT} ]; then
  Analyze \
    ${EDEP_OUT} \
    ${RECO_OUT}
fi
```

The edep-sim macro template is:

```
################################
# Add the GENIE events
################################

## Replace this with the name of a GENIE rooTracker file
/generator/kinematics/rooTracker/input __INPUT__

## Use the T2K rooTracker input format.  This is directly supported by GENIE.
/generator/kinematics/set rooTracker

## Distribute the events based on the density of the material.  When done
##   this way, the composition of the detector is ignored, so it's not
##   a good way for physics, but it's OK for an example since you don't
##   need to syncronize the GENIE and EDEPSIM geometries.
#/generator/position/density/sample DetEnclosure_lv
#/generator/position/set density

#/generator/count/fixed/number 1000
#/generator/count/set fixed

#/generator/time/spill/start 0.0 ns
#/generator/time/spill/bunchCount 1000
#/generator/time/spill/bunchSep 10.0 ns
#/generator/time/spill/bunchLength 5.0 ns

## Make sure EDEPSIM updates the kinematics generator.
/generator/add

/run/initialize
/run/beamOn __NEV__
```

</details>

# Data format

The description of the data format can be found [here](https://github.com/DUNE-ND-SAND/sand-stt/wiki/Data-Model)


# Contribute

- The code format can be find [here](https://github.com/DUNE-ND-SAND/sand-stt/wiki/Code-Formatting)
- The developing scheme is described [HowToDevelop.pdf](https://baltig.infn.it/dune/sand-reco/-/wikis/uploads/8b897fb0ea753ef767b96312bdf9ccac/HowToDevelop.pdf)

# Support

For any communication, please refer to [DUNE-ND-SAND-SOFTWARE@fnal.gov](mailto:DUNE-ND-SAND-SOFTWARE@fnal.gov)
