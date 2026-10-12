# Repeating experiments

## Setup

In the main CMR repository, run
```
mkdir build-release
cd build-release
cmake .. -DGENERATORS=1
```
to compile CMR with instance generators.
The executables `cmr-tu`, `cmr-generate-network`, etc. should then reside in the current (build) directory.

### Install predecessor code

In a separate location, one also needs to install the unimodularity test (latest version 1.2h) from http://matthiaswalter.org/TUtest/:
```
wget http://matthiaswalter.org/TUtest/unimodularity-library-1.2h.tar.gz
tar xzf unimodularity-library-1.2h.tar.gz
cd unimodularity-library-1.2h/
./configure
make
```
The executable `unimodularity-test` should then reside in `src/`.

## Experiments

### Configuration

Before running the experiments, one has to centrally set up the directories.
The file `local.py` should set

 - `BUILD_DIRECTORY` to the full path of the `build-release` directory created above,
 - `TEMP_STORAGE` to the full path of temporary instance storage (e.g., if a cluster is used that has local storage).
 - `RESULTS_STORAGE` to the directory containing the results.
 - `UNIMOD_DIRECTORY` to the ful path of the `src` folder in the unimodularity-library.

Except for the latter, reasonable default values are given.

### Infrastructure

The experiments are separated into so-called /runs/, which can be individually executed in any order besides some dependencies.
A run is specified as `CLASS-SIZE-SAMPLE-ALGO` with these meanings:

 - `CLASS` refers to a certain class of matrices such as `modwheel`, `network`, `rndcamion` or `identity`.
 - `SIZE` is an indicator of the matrix size, e.g., `00019x00019` for a 19-by-19 matrix.
 - `SAMPLE` indicates the number of the sample for experiments that a repeated a number of times.
 - `ALGO` refers to an algorithm of the form `cmr*` or `unimod*`, e.g., `cmrcert` for CMR's decomposition algorithm with additional search for a NO-certificate (non-TU submatrix). Moreover, the value `gen` refers to instance generation.

The set of runs that are configured can be listed by executing `python3 runs.py`.
In order to modify the runs or options like the time limit, `runs.py` can be modified.

Once the runs are determined, one can carry out experiments by means of the `execute.py` file.
It contains the actual logic, i.e., the calls of the `cmr-*` tools, etc., and logic for parallelization.
One can run multiple instances at once, e.g., by calling `python3 execute.py` multiple times.
The processes coordinate their work using files like `results/RUN.lock` and `results/RUN.progress`, where `RUN` refers to the respective run.
If one run is being processed, `execute.py` simply picks another one.
For technical reasons the lock file is not deleted.
The progress file contains information about whether the run is still being processed or complete.
By executing `python3 runs.py '.*' status` the status of all runs can be inspected.
The regular expression parameter can be adapted to only report about a subset of runs.

The file `execute.sbatch` contains an example setup for SLURM-based configuration with 10 parallel processes.

### Evaluation

Statistics about the results can be gathered by
```
python3 eval_modwheel.py  # Difficult non-TU matrices
python3 eval_network.py   # Network matrices
python3 eval_rndcamion.py # Random 0/1-matrices that were signed using the cmr-camion tool.
```
all of which accept an optional regular expression parameter for reporting only about a subset of runs.
If that is followed by `boxplot` then LaTeX code for an actual box plot is printed instead.

