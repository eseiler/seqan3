// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/utility/range/to.hpp>
#include <seqan3/utility/views/join_with.hpp>

struct options
{
    size_t const sequence_length;
    bool const has_repeats;
    size_t const number_of_reads;
    size_t const read_length;
    double const prob_insertion;
    double const prob_deletion;
    uint8_t const simulated_errors;
    uint8_t const searched_errors;
    uint8_t const strata;
    double const stddev{0};
    uint32_t repeats{20};
};

template <seqan3::alphabet alphabet_t>
void mutate_substitution(std::vector<alphabet_t> & seq, size_t const pos, uint8_t alphabet_rank)
{
    alphabet_t & cbase = seq[pos];
    if (alphabet_rank >= seqan3::to_rank(cbase))
        ++alphabet_rank;
    cbase.assign_rank(alphabet_rank);
}

template <seqan3::alphabet alphabet_t>
void mutate_insertion(std::vector<alphabet_t> & seq, size_t const pos, uint8_t const alphabet_rank)
{
    seq.insert(std::ranges::begin(seq) + pos, alphabet_t{}.assign_rank(alphabet_rank));
}

template <seqan3::alphabet alphabet_t>
void mutate_deletion(std::vector<alphabet_t> & seq, size_t const pos)
{
    seq.erase(std::ranges::begin(seq) + pos);
}

template <seqan3::alphabet alphabet_t>
std::vector<std::vector<alphabet_t>> generate_reads(std::vector<alphabet_t> const & ref,
                                                    size_t const number_of_reads,
                                                    size_t const read_length,
                                                    uint8_t simulated_errors,
                                                    double const prob_insertion,
                                                    double const prob_deletion,
                                                    double const stddev = 0,
                                                    size_t const seed = 0)
{
    std::vector<std::vector<alphabet_t>> reads;
    std::mt19937_64 gen{seed};

    auto error_count_distribution = [simulated_errors, stddev](auto & generator) -> uint8_t
    {
        if (stddev <= 0)
            return simulated_errors;

        std::normal_distribution<> distribution{static_cast<double>(simulated_errors), stddev};
        return std::abs(std::rint(distribution(generator)));
    };

    // mutation distributions
    std::uniform_real_distribution<double> mutation_type_prob{0.0, 1.0};
    // position
    std::uniform_int_distribution<size_t> random_mutation_pos{0, read_length - 1};
    // substitution
    std::uniform_int_distribution<uint8_t> dis_alpha_short{0, seqan3::alphabet_size<alphabet_t> - 2};
    // insertion
    std::uniform_int_distribution<uint8_t> dis_alpha{0, seqan3::alphabet_size<alphabet_t> - 1};

    for (size_t i = 0; i < number_of_reads; ++i)
    {
        // simulate concrete error number or use normal distribution
        simulated_errors = error_count_distribution(gen);

        std::uniform_int_distribution<size_t> random_read_pos{0,
                                                              std::ranges::size(ref) - read_length - simulated_errors};
        size_t rpos = random_read_pos(gen);
        std::vector<alphabet_t> read_tmp{std::ranges::begin(ref) + rpos,
                                         std::ranges::begin(ref) + rpos + read_length + simulated_errors};

        // generate simulated_errors many unique random mutation positions
        std::set<size_t> mutation_positions;
        if (read_length > simulated_errors)
        {
            while (mutation_positions.size() < simulated_errors)
                mutation_positions.insert(random_mutation_pos(gen));
        }
        else
        {
            for (size_t i = 0; i < simulated_errors; ++i)
                mutation_positions.insert(i);
        }

        for (std::set<size_t>::iterator pos_it = mutation_positions.begin(); pos_it != mutation_positions.end();
             ++pos_it)
        {
            size_t ppos = *pos_it;
            double prob = mutation_type_prob(gen);
            // Substitution
            if (prob_insertion + prob_deletion < prob)
                mutate_substitution(read_tmp, ppos, dis_alpha_short(gen));
            // Insertion
            else if (prob_insertion < prob)
                mutate_insertion(read_tmp, ppos, dis_alpha(gen));
            // Deletion
            else
                mutate_deletion(read_tmp, ppos);
        }

        read_tmp.erase(std::ranges::begin(read_tmp) + read_length, std::ranges::end(read_tmp));
        reads.push_back(read_tmp);
    }

    return reads;
}

template <typename alphabet_t>
std::vector<alphabet_t> generate_repeating_sequence(size_t const template_length = 5000,
                                                    size_t const repeats = 20,
                                                    double const template_fraction = 1,
                                                    size_t const seed = 0)
{
    std::vector<alphabet_t> seq_template = seqan3::test::generate_sequence<alphabet_t>(template_length, 0, seed);

    // copy substrings of length len from seq_template mutate and concatenate them
    size_t len = std::round(template_length * template_fraction);
    uint8_t simulated_errors = 5;
    len = (len + simulated_errors > template_length) ? template_length - simulated_errors : len;

    return generate_reads(seq_template, repeats, len, simulated_errors, 0.15, 0.15) | std::views::all | std::views::join
         | seqan3::ranges::to<std::vector>();
}

//============================================================================
//  undirectional; trivial_search, collection, dna4, all-mapping
//============================================================================

// Note: We force the actual computation of the search() by going through the lazy algorithm result range (input_range)
//       using std::ranges::distance within the for loop.

void unidirectional_search_all_collection(benchmark::State & state, options && o)
{
    size_t set_size = 10;
    std::vector<std::vector<seqan3::dna4>> collection;
    std::vector<std::vector<seqan3::dna4>> reads;
    for (size_t i = 0; i < set_size; ++i)
    {
        collection.push_back(seqan3::test::generate_sequence<seqan3::dna4>(o.sequence_length, 0, i));
        std::vector<std::vector<seqan3::dna4>> seq_reads = generate_reads(collection.back(),
                                                                          o.number_of_reads,
                                                                          o.read_length,
                                                                          o.simulated_errors,
                                                                          o.prob_insertion,
                                                                          o.prob_deletion,
                                                                          o.stddev,
                                                                          i);
        std::ranges::move(seq_reads, std::back_inserter(reads));
    }

    seqan3::fm_index index{collection};
    seqan3::configuration cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{o.searched_errors}};

    size_t sum{};
    for (auto _ : state)
    {
        auto results = search(reads, index, cfg);
        sum += std::ranges::distance(results);
    }
    benchmark::DoNotOptimize(sum);
}

//============================================================================
//  undirectional; trivial_search, single, dna4, all-mapping
//============================================================================

void unidirectional_search_all(benchmark::State & state, options && o)
{
    std::vector<seqan3::dna4> ref =
        (o.has_repeats)
            ? generate_repeating_sequence<seqan3::dna4>(2 * o.sequence_length / o.repeats, o.repeats, 0.5, 0)
            : seqan3::test::generate_sequence<seqan3::dna4>(o.sequence_length, 0, 0);

    seqan3::fm_index index{ref};
    std::vector<std::vector<seqan3::dna4>> reads = generate_reads(ref,
                                                                  o.number_of_reads,
                                                                  o.read_length,
                                                                  o.simulated_errors,
                                                                  o.prob_insertion,
                                                                  o.prob_deletion,
                                                                  o.stddev);
    seqan3::configuration cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{o.searched_errors}};

    size_t sum{};
    for (auto _ : state)
    {
        auto results = search(reads, index, cfg);
        sum += std::ranges::distance(results);
    }
    benchmark::DoNotOptimize(sum);
}

//============================================================================
//  bidirectional; trivial_search, single, dna4, all-mapping
//============================================================================

void bidirectional_search_all(benchmark::State & state, options && o)
{
    std::vector<seqan3::dna4> ref =
        (o.has_repeats)
            ? generate_repeating_sequence<seqan3::dna4>(2 * o.sequence_length / o.repeats, o.repeats, 0.5, 0)
            : seqan3::test::generate_sequence<seqan3::dna4>(o.sequence_length, 0, 0);

    seqan3::bi_fm_index index{ref};
    std::vector<std::vector<seqan3::dna4>> reads = generate_reads(ref,
                                                                  o.number_of_reads,
                                                                  o.read_length,
                                                                  o.simulated_errors,
                                                                  o.prob_insertion,
                                                                  o.prob_deletion,
                                                                  o.stddev);
    seqan3::configuration cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{o.searched_errors}};

    size_t sum{};
    for (auto _ : state)
    {
        auto results = search(reads, index, cfg);
        sum += std::ranges::distance(results);
    }
    benchmark::DoNotOptimize(sum);
}

//============================================================================
//  undirectional; trivial_search, single, dna4, stratified-all-mapping
//============================================================================

void unidirectional_search_stratified(benchmark::State & state, options && o)
{
    std::vector<seqan3::dna4> ref =
        (o.has_repeats)
            ? generate_repeating_sequence<seqan3::dna4>(2 * o.sequence_length / o.repeats, o.repeats, 0.5, 0)
            : seqan3::test::generate_sequence<seqan3::dna4>(o.sequence_length, 0, 0);

    seqan3::fm_index index{ref};
    std::vector<std::vector<seqan3::dna4>> reads = generate_reads(ref,
                                                                  o.number_of_reads,
                                                                  o.read_length,
                                                                  o.simulated_errors,
                                                                  o.prob_insertion,
                                                                  o.prob_deletion,
                                                                  o.stddev);
    seqan3::configuration cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{o.searched_errors}}
                              | seqan3::search_cfg::hit_strata{o.strata};

    size_t sum{};
    for (auto _ : state)
    {
        auto results = search(reads, index, cfg);
        sum += std::ranges::distance(results);
    }
    benchmark::DoNotOptimize(sum);
}

//============================================================================
//  bidirectional; trivial_search, single, dna4, stratified-all-mapping
//============================================================================

void bidirectional_search_stratified(benchmark::State & state, options && o)
{
    std::vector<seqan3::dna4> ref =
        (o.has_repeats)
            ? generate_repeating_sequence<seqan3::dna4>(2 * o.sequence_length / o.repeats, o.repeats, 0.5, 0)
            : seqan3::test::generate_sequence<seqan3::dna4>(o.sequence_length, 0, 0);

    seqan3::bi_fm_index index{ref};
    std::vector<std::vector<seqan3::dna4>> reads = generate_reads(ref,
                                                                  o.number_of_reads,
                                                                  o.read_length,
                                                                  o.simulated_errors,
                                                                  o.prob_insertion,
                                                                  o.prob_deletion,
                                                                  o.stddev);
    seqan3::configuration cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{o.searched_errors}}
                              | seqan3::search_cfg::hit_strata{o.strata};

    size_t sum{};
    for (auto _ : state)
    {
        auto results = search(reads, index, cfg);
        sum += std::ranges::distance(results);
    }
    benchmark::DoNotOptimize(sum);
}

#ifndef NDEBUG
inline constexpr size_t small_size = 1000;
inline constexpr size_t medium_size = 5000;
inline constexpr size_t big_size = 10000;
#else
inline constexpr size_t small_size = 10000;
inline constexpr size_t medium_size = 50000;
inline constexpr size_t big_size = 100000;
#endif // NDEBUG

BENCHMARK_CAPTURE(unidirectional_search_all_collection,
                  highErrorReadsSearch0,
                  options{small_size, false, 10, 50, 0.18, 0.18, 0, 0, 0, 1.75});
BENCHMARK_CAPTURE(unidirectional_search_all_collection,
                  highErrorReadsSearch1,
                  options{small_size, false, 10, 50, 0.18, 0.18, 0, 1, 0, 1.75});
BENCHMARK_CAPTURE(unidirectional_search_all_collection,
                  highErrorReadsSearch2,
                  options{small_size, false, 10, 50, 0.18, 0.18, 0, 2, 0, 1.75});
BENCHMARK_CAPTURE(unidirectional_search_all_collection,
                  highErrorReadsSearch3,
                  options{small_size, false, 10, 50, 0.18, 0.18, 0, 3, 0, 1.75});

BENCHMARK_CAPTURE(unidirectional_search_all,
                  lowErrorReadsSearch3,
                  options{big_size, false, 50, 50, 0.18, 0.18, 0, 3, 0, 1});
BENCHMARK_CAPTURE(unidirectional_search_all,
                  highErrorReadsSearch0,
                  options{big_size, false, 50, 50, 0.18, 0.18, 0, 0, 0, 1.75});
BENCHMARK_CAPTURE(unidirectional_search_all,
                  highErrorReadsSearch1,
                  options{big_size, false, 50, 50, 0.18, 0.18, 0, 1, 1, 1.75});
BENCHMARK_CAPTURE(unidirectional_search_all,
                  highErrorReadsSearch2,
                  options{big_size, false, 50, 50, 0.18, 0.18, 0, 2, 2, 1.75});
BENCHMARK_CAPTURE(unidirectional_search_all,
                  highErrorReadsSearch3,
                  options{big_size, false, 50, 50, 0.18, 0.18, 0, 3, 3, 1.75});
BENCHMARK_CAPTURE(unidirectional_search_all,
                  highErrorReadsSearch0Rep,
                  options{big_size, true, 50, 50, 0.18, 0.18, 0, 0, 0, 1.75});
BENCHMARK_CAPTURE(unidirectional_search_all,
                  highErrorReadsSearch1Rep,
                  options{big_size, true, 50, 50, 0.18, 0.18, 0, 1, 1, 1.75});
BENCHMARK_CAPTURE(unidirectional_search_all,
                  highErrorReadsSearch2Rep,
                  options{big_size, true, 50, 50, 0.18, 0.18, 0, 2, 2, 1.75});
BENCHMARK_CAPTURE(unidirectional_search_all,
                  highErrorReadsSearch3Rep,
                  options{big_size, true, 50, 50, 0.18, 0.18, 0, 3, 3, 1.75});
BENCHMARK_CAPTURE(unidirectional_search_all,
                  highErrorReadsSearch3Rep,
                  options{big_size, true, 50, 50, 0.30, 0.30, 0, 3, 3, 1.75});

BENCHMARK_CAPTURE(bidirectional_search_all,
                  lowErrorReadsSearch3,
                  options{big_size, false, 50, 50, 0.18, 0.18, 0, 3, 0, 1});
BENCHMARK_CAPTURE(bidirectional_search_all,
                  highErrorReadsSearch0,
                  options{big_size, false, 50, 50, 0.18, 0.18, 0, 0, 0, 1.75});
BENCHMARK_CAPTURE(bidirectional_search_all,
                  highErrorReadsSearch1,
                  options{big_size, false, 50, 50, 0.18, 0.18, 0, 1, 1, 1.75});
BENCHMARK_CAPTURE(bidirectional_search_all,
                  highErrorReadsSearch2,
                  options{big_size, false, 50, 50, 0.18, 0.18, 0, 2, 2, 1.75});
BENCHMARK_CAPTURE(bidirectional_search_all,
                  highErrorReadsSearch3,
                  options{big_size, false, 50, 50, 0.18, 0.18, 0, 3, 3, 1.75});
BENCHMARK_CAPTURE(bidirectional_search_all,
                  highErrorReadsSearch0Rep,
                  options{big_size, true, 50, 50, 0.18, 0.18, 0, 0, 0, 1.75});
BENCHMARK_CAPTURE(bidirectional_search_all,
                  highErrorReadsSearch1Rep,
                  options{big_size, true, 50, 50, 0.18, 0.18, 0, 1, 1, 1.75});
BENCHMARK_CAPTURE(bidirectional_search_all,
                  highErrorReadsSearch2Rep,
                  options{big_size, true, 50, 50, 0.18, 0.18, 0, 2, 2, 1.75});
BENCHMARK_CAPTURE(bidirectional_search_all,
                  highErrorReadsSearch3Rep,
                  options{big_size, true, 50, 50, 0.18, 0.18, 0, 3, 3, 1.75});
BENCHMARK_CAPTURE(bidirectional_search_all,
                  highErrorReadsSearch3Rep,
                  options{big_size, true, 50, 50, 0.30, 0.30, 0, 3, 3, 1.75});

BENCHMARK_CAPTURE(unidirectional_search_stratified,
                  lowErrorReadsSearch3Strata0Rep,
                  options{medium_size, true, 50, 50, 0.18, 0.18, 0, 3, 0, 1});
BENCHMARK_CAPTURE(unidirectional_search_stratified,
                  lowErrorReadsSearch3Strata1Rep,
                  options{medium_size, true, 50, 50, 0.18, 0.18, 0, 3, 1, 1});
BENCHMARK_CAPTURE(unidirectional_search_stratified,
                  lowErrorReadsSearch3Strata2Rep,
                  options{medium_size, true, 50, 50, 0.18, 0.18, 0, 3, 2, 1});
BENCHMARK_CAPTURE(unidirectional_search_stratified,
                  highErrorReadsSearch3Strata0Rep,
                  options{medium_size, true, 50, 50, 0.30, 0.30, 0, 3, 0, 1.75});
BENCHMARK_CAPTURE(unidirectional_search_stratified,
                  highErrorReadsSearch3Strata1Rep,
                  options{medium_size, true, 50, 50, 0.30, 0.30, 0, 3, 1, 1.75});
BENCHMARK_CAPTURE(unidirectional_search_stratified,
                  highErrorReadsSearch3Strata2Rep,
                  options{medium_size, true, 50, 50, 0.30, 0.30, 0, 3, 2, 1.75});
BENCHMARK_CAPTURE(unidirectional_search_stratified,
                  highErrorReadsSearch3Strata2RepLong,
                  options{big_size, true, 50, 50, 0.30, 0.30, 0, 3, 2, 1.75});

BENCHMARK_CAPTURE(bidirectional_search_stratified,
                  lowErrorReadsSearch3Strata0Rep,
                  options{medium_size, true, 50, 50, 0.18, 0.18, 0, 3, 0, 1});
BENCHMARK_CAPTURE(bidirectional_search_stratified,
                  lowErrorReadsSearch3Strata1Rep,
                  options{medium_size, true, 50, 50, 0.18, 0.18, 0, 3, 1, 1});
BENCHMARK_CAPTURE(bidirectional_search_stratified,
                  lowErrorReadsSearch3Strata2Rep,
                  options{medium_size, true, 50, 50, 0.18, 0.18, 0, 3, 2, 1});
BENCHMARK_CAPTURE(bidirectional_search_stratified,
                  highErrorReadsSearch3Strata0Rep,
                  options{medium_size, true, 50, 50, 0.30, 0.30, 0, 3, 0, 1.75});
BENCHMARK_CAPTURE(bidirectional_search_stratified,
                  highErrorReadsSearch3Strata1Rep,
                  options{medium_size, true, 50, 50, 0.30, 0.30, 0, 3, 1, 1.75});
BENCHMARK_CAPTURE(bidirectional_search_stratified,
                  highErrorReadsSearch3Strata2Rep,
                  options{medium_size, true, 50, 50, 0.30, 0.30, 0, 3, 2, 1.75});
BENCHMARK_CAPTURE(bidirectional_search_stratified,
                  highErrorReadsSearch3Strata2RepLong,
                  options{big_size, true, 50, 50, 0.30, 0.30, 0, 3, 2, 1.75});

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
