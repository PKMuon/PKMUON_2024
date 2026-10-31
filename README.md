# PKMuon Geant4 Simulation Program 2024

RPC response reference tools and the measured-resolution coordinate baseline are
documented in [rpc_response/README.md](rpc_response/README.md). See
[validation results](rpc_response/VALIDATION.md) for the first Ubuntu run.
中文[程序结构与物理原理说明](rpc_response/PROGRAM_STRUCTURE_AND_PHYSICS.md)包含各模块、数据接口、物理模型和适用范围。

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
