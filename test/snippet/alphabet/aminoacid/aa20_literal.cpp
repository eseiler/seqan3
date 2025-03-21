// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/aminoacid/aa20.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aa20_vector sequence1{"ACGTTA"_aa20};
    seqan3::aa20_vector sequence2 = "ACGTTA"_aa20;
    auto sequence3 = "ACGTTA"_aa20;
}
