# PKMuon Geant4 Simulation Program 2024

Contact \<seeson@pku.edu.cn\> on any problem.

Setup ROOT and Geant4 environment before everything.

Compile at the first time:

```bash
cd PKMUON_2024
mkdir build
cd build
cmake ..
make -j8
```

Recompile (needed on any modification on source files (`.cc`)) afterwards:

```bash
cd PKMUON_2024/build
make -j8
```

Launch GUI:

```bash
cd PKMUON_2024/build
./muPos
```

Launch CryMu simulation (with multi-process parallelism):

```bash
cd PKMUON_2024/build
./run_cry.sh
```

Analyze data from the previous step:

```bash
cd PKMUON_2024/analysis
./analysis.sh
```
