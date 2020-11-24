Build
- make sure you have *root-config* in your ${PATH} environment variable
- open Makefile and customize the location of *Edep-sim* installation
- Then:

```
$ make
```

Before running application
- To have dictionaties of the structs loaded at run time:
```
$ source scripts/env.sh
```

Using with ROOT
- To load dictionaries of structs:
```
root [0] gSystem->Load("libStruct.so")
```

Digitization
- Create digits of STT and cells of calorimeter

```
$ Digitize <input file> <output file>
```

Reconstruction
- Track find and fit of STT track
- Clustering of calorimeter cells

```
$ Reconstruct <input file>
```

Analysis
- Evaluate parameters of particles
- Evaluate neutrino energy

```
Analyze <input file>
```

The description of the data format can be found [here](https://github.com/DUNE-ND-SAND/sand-stt/wiki/Data-Model)

The code format can be find [here](https://github.com/DUNE-ND-SAND/sand-stt/wiki/Code-Formatting)
