// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@.cpp.in

//![main]
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::rna15 letter1 = 'C'_dna15; // implicitly converted
    seqan3::rna15 letter2{};
    letter2 = 'C'_dna15; // implicitly converted
}
//![main]
