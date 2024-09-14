# Description

A program for sequential and parallel (OMP) calculation of DOS (probability density of states) for a system of dipoles. There is a consideration of periodic boundary conditions.

# How to build

```
git clone --recurse-submodules https://github.com/uyras/bruteforce.git
cd bruteforce
mkdir build && cd build
cmake ..
cmake --build .
```

# How to use

`./bruteForce_seq` - is the sequential code. Run it with `--help` to see the command-line parameters.

Run `./bruteForce_omp --help` for parallel code. You can manage the number of threads by changing env. variable `OMP_NUM_THREADS`.

Example: `OMP_NUM_THREADS=5 ./bruteForce_omp`. By default it utilise all avaliable CPUs.
