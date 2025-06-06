// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_views.cpp.in

//![main]
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/utility/views/convert.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::rna5_vector vector = "ACG"_rna5;

    auto dna5_view = vector | seqan3::views::convert<seqan3::dna5>;

    for (auto && chr : dna5_view) // converts lazily on-the-fly
    {
        static_assert(std::same_as<decltype(chr), seqan3::dna5 &&>);
    }
}
//![main]
