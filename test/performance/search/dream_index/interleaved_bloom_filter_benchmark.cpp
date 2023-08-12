// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#define OLD__ 0

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/utility/range/to.hpp>
#include <seqan3/utility/views/zip.hpp>

inline benchmark::Counter hashes_per_second(size_t const count)
{
    return benchmark::Counter(count, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1000);
}

#if 1
static void arguments(benchmark::internal::Benchmark * b)
{
    // Total size: 1MiB
    // bins, bin_size, hash_num, sequence_length
    b->Args({64, 1LL << 17, 2, 1LL << 17});
    b->Args({128, 1LL << 16, 2, 1LL << 17});
    b->Args({192, 1LL << 16, 2, 1LL << 17});
    b->Args({256, 1LL << 15, 2, 1LL << 17});
    b->Args({1024, 1LL << 10, 2, 1LL << 17});
}
#else
static void arguments(benchmark::internal::Benchmark * b)
{
    // Total size: 1GiB
    // bins, bin_size, hash_num, sequence_length
    b->Args({64, 1LL << 27, 2, 1LL << 27});
    b->Args({8192, 1LL << 20, 2, 1LL << 27});
}
#endif

template <typename ibf_type>
auto set_up(size_t bins, size_t bits, size_t hash_num, size_t sequence_length)
{
    auto bin_indices = seqan3::test::generate_numeric_sequence<size_t>(sequence_length, 0u, bins - 1);
    auto hash_values = seqan3::test::generate_numeric_sequence<size_t>(sequence_length);
    seqan3::interleaved_bloom_filter tmp_ibf(seqan3::bin_count{bins},
                                             seqan3::bin_size{bits},
                                             seqan3::hash_function_count{hash_num});

    ibf_type ibf{std::move(tmp_ibf)};

    return std::make_tuple(bin_indices, hash_values, ibf);
}

template <typename ibf_type>
void bulk_contains_benchmark(::benchmark::State & state)
{
    auto && [bin_indices, hash_values, ibf] =
        set_up<ibf_type>(state.range(0), state.range(1), state.range(2), state.range(3));

    for (auto [hash, bin] : seqan3::views::zip(hash_values, bin_indices))
        ibf.emplace(hash, seqan3::bin_index{bin});

    auto agent = ibf.membership_agent();
    for (auto _ : state)
    {
        for (auto hash : hash_values)
        {
            [[maybe_unused]] auto & res = agent.bulk_contains(hash);
            benchmark::ClobberMemory();
        }
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
}

// template <typename ibf_type>
// void bulk_count_benchmark(::benchmark::State & state)
// {
//     auto && [bin_indices, hash_values, ibf] =
//         set_up<ibf_type>(state.range(0), state.range(1), state.range(2), state.range(3));
//     (void)bin_indices;

//     auto agent = ibf.counting_agent();
//     for (auto _ : state)
//     {
//         [[maybe_unused]] auto & res = agent.bulk_count(hash_values);
//     }

//     state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
// }

BENCHMARK_TEMPLATE(bulk_contains_benchmark, seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>)
    ->Apply(arguments);

// BENCHMARK_TEMPLATE(bulk_count_benchmark, seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>)
//     ->Apply(arguments);

BENCHMARK_MAIN();
