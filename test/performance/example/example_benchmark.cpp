// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstring>
#include <vector>

#include <seqan3/test/performance/units.hpp>

static void vector_copy_benchmark(benchmark::State & state)
{
    std::vector<int> x = {15, 13, 12, 10};
    for (auto _ : state)
        std::vector<int> copy{x};
}

static void memcpy_benchmark(benchmark::State & state)
{
    size_t size = state.range(0);
    std::vector<char> src(size, '-');
    std::vector<char> dst(size);

    for (auto _ : state)
        std::copy(src.begin(), src.end(), dst.begin());

    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(size);
}

static void copy_benchmark(benchmark::State & state)
{
    unsigned size = state.range(0);
    char * src = new char[size];
    char * dst = new char[size];

    memset(src, '-', size);
    for (auto _ : state)
        std::copy_n(src, size, dst);

    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(size);
    delete[] src;
    delete[] dst;
}

// Register the function as a benchmark
BENCHMARK(vector_copy_benchmark);

BENCHMARK(memcpy_benchmark)->Arg(8)->Arg(64)->Arg(512);
BENCHMARK(memcpy_benchmark)->Range(4, 4 << 5);
BENCHMARK(copy_benchmark)->RangeMultiplier(2)->Range(4, 4 << 5);

BENCHMARK_MAIN();
