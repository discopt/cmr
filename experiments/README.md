# Repeating experiments

## Setup

In the main repository, run
```
mkdir build-release
cd build-release
cmake .. -DGENERATORS=1
```
to compile CMR with instance generators.

## Experiments

Before running the experiments, one has to centrally set up the directories.
The file `config.py` should set
- `BUILD_DIRECTORY` to the full path of the `build-release` directory created above,
- `UNIMOD_DIRECTORY` to the ful path of the `src` folder in the predecessor unimodularity-library.
- `LOCAL_STORAGE` to the full path of temporary instance storage.
- `ALGORITHMS` does not need to be modified.

### Difficult non-totally unimodular matrices

First, create the results directory.
```
mkdir modwheel
```
Second
```
python3 gen_modwheel.py
```
generates the instances.
Alternatively, `gen_modwheel.sbatch` can be used to schedule it as a job in a slurm-controlled cluster environment.
Third,
```
python3 run_modwheel.py <SEED>
```
should be run for different seed values `1`, `2`, ..., `10`.
Alternatively, `run_modwheel.sbatch` can be used to schedule it as a job in a slurm-controlled cluster environment.
