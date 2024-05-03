// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_and_window_hash.hpp>
#include <seqan3/search/views/minimiser_hash_and_positions.hpp>

using seqan3::operator""_dna4;

static seqan3::dna4_vector const text{"TCATCAGTAGCTACAATACG"_dna4};
static constexpr size_t const minimiser_size = 2;
static constexpr size_t const window_size = 5;

seqan3::dna4_vector kmer_to_string(uint64_t kmer, size_t const kmer_size)
{
    seqan3::dna4_vector result(kmer_size);
    for (size_t i = 0; i < kmer_size; ++i)
    {
        result[kmer_size - 1 - i].assign_rank(kmer & 0b11);
        kmer >>= 2;
    }
    return result;
}

void minimiser_hash_and_positions()
{
    auto view =
        bsc::views::minimiser_hash_and_positions({.minimiser_size = minimiser_size, .window_size = window_size});

    seqan3::debug_stream << "minimiser,position,occurrences\n";

    for (auto && result : text | view)
    {
        seqan3::debug_stream << kmer_to_string(result.minimiser_value, minimiser_size) << ',' << result.range_position
                             << ',' << result.occurrences << '\n';
    }

    seqan3::debug_stream << '\n';

    // AT,0,2
    // AG,2,7
    // AC,9,2
    // AA,11,4
    // AC,15,1
}

void minimiser_and_window_hash()
{
    auto view = bsc::views::minimiser_and_window_hash({.minimiser_size = minimiser_size, .window_size = window_size});

    seqan3::debug_stream << "minimiser,window\n";

    for (auto && result : text | view)
    {
        seqan3::debug_stream << kmer_to_string(result.minimiser_value, minimiser_size) << ','
                             << kmer_to_string(result.window_value, window_size) << '\n';
    }

    seqan3::debug_stream << '\n';

    // AT,TCATC
    // AG,ATCAG
    // AG,GTAGC
    // AC,GCTAC
    // AA,TACAA
    // AC,ATACG
}

int main()
{
    seqan3::debug_stream << "Text\n" << text << "\n01234567890123456789\n\n";
    minimiser_hash_and_positions();
    minimiser_and_window_hash();
}
