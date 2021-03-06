#!/bin/bash

export testfile=performance_testing
export inputfile=performance_testing_input
export branch1=develop
export branch2=optimize_inc_beta

clang-format-5.0 -i $testfile.cpp
clang-format-5.0 -i $inputfile.cpp

benchmark_branch () {
    local branch=$1
    cd ../math
    git checkout --quiet $branch
    git rev-parse --abbrev-ref HEAD
    cd ..
    g++ math-microbenchmark/$testfile.cpp math-microbenchmark/$inputfile.cpp -DBRANCH=$branch -std=c++1y -D_REENTRANT -Wno-sign-compare -isystem benchmark/include -Lbenchmark/build/src -lbenchmark -lpthread -I math/lib/tbb_2019_U8/include -O3 -I math/ -I math/lib/eigen_3.3.3 -I math/lib/boost_1.72.0 -I math/lib/sundials_4.1.0/include -DBOOST_DISABLE_ASSERTS  -Wl,-L,"math/lib/tbb" -Wl,-rpath,"math/lib/tbb" math/lib/tbb/libtbb.so.2  math/lib/sundials_4.1.0/lib/libsundials_nvecserial.a math/lib/sundials_4.1.0/lib/libsundials_cvodes.a math/lib/sundials_4.1.0/lib/libsundials_idas.a math/lib/sundials_4.1.0/lib/libsundials_kinsol.a -o math-microbenchmark/$testfile
    chmod +777 math-microbenchmark/$testfile
    math-microbenchmark/$testfile --benchmark_report_aggregates_only
    rm math-microbenchmark/$testfile
    cd math-microbenchmark
}

benchmark_branch $branch1
benchmark_branch $branch2
