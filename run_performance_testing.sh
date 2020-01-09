#!/bin/bash

export testfile=performance_testing
export inputfile=performance_testing_input
export branch1=develop
export branch2=optimize_incomplete_beta
export branch3=optimize_incomplete_beta~
export pathtotop=../

clang-format-5.0 -i $testfile.cpp
clang-format-5.0 -i $inputfile.cpp

benchmark_branch () {
    local branch=$1
    cd $pathtotop
    cd math
    git checkout --quiet $branch
    git rev-parse --abbrev-ref HEAD
    cd $pathtotop
    g++ math-microbenchmark/$testfile.cpp math-microbenchmark/$inputfile.cpp -DBRANCH=$branch -std=c++1y -D_REENTRANT -Wno-sign-compare -isystem benchmark/include -Lbenchmark/build/src -lbenchmark -lpthread -I math/lib/tbb_2019_U8/include -O3 -I math/ -I math/lib/eigen_3.3.3 -I math/lib/boost_1.72.0 -I math/lib/sundials_4.1.0/include -DBOOST_DISABLE_ASSERTS  -Wl,-L,"math/lib/tbb" -Wl,-rpath,"math/lib/tbb" math/lib/tbb/libtbb.so.2  math/lib/sundials_4.1.0/lib/libsundials_nvecserial.a math/lib/sundials_4.1.0/lib/libsundials_cvodes.a math/lib/sundials_4.1.0/lib/libsundials_idas.a math/lib/sundials_4.1.0/lib/libsundials_kinsol.a -o $testfile
    chmod +777 $testfile
    ./$testfile --benchmark_report_aggregates_only
    rm $testfile
}

benchmark_branch $branch1
benchmark_branch $branch2
benchmark_branch $branch3
