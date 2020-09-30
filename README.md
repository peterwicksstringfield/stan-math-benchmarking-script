Don't use this. See [stan-dev/perf-math](https://github.com/stan-dev/perf-math).

Kludgy thing that benchmarks some stuff in the math library.

Run like this:
./run_performance_testing.sh

Requires [benchmark](https://github.com/google/benchmark).

```
/path/to/top/
    math/
        .git (branch1 and branch2)
    stan-math-benchmarking-script/
        .git
        testfile
        inputfile
        run_performance_testing.sh
    benchmark/
        .git
```

Benchmarks in testfile. Input in inputfile. Separate translation unit to
prevent spurious constant propagation.

Choose the branches to compare in the shell script.
