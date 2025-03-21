// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/core/debug_stream/range.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "nucleotide_test_template.hpp"

using seqan3::operator""_dna15;

INSTANTIATE_TYPED_TEST_SUITE_P(dna15, alphabet, seqan3::dna15, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna15, semi_alphabet_test, seqan3::dna15, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna15, alphabet_constexpr, seqan3::dna15, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna15, semi_alphabet_constexpr, seqan3::dna15, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna15, nucleotide, seqan3::dna15, );

TEST(dna15, to_char_assign_char)
{
    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('A')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('C')), 'C');
    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('G')), 'G');

    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('U')), 'T');
    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('T')), 'T');

    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('R')), 'R');
    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('Y')), 'Y');
    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('S')), 'S');
    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('W')), 'W');
    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('K')), 'K');
    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('M')), 'M');
    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('B')), 'B');
    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('D')), 'D');
    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('H')), 'H');
    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('V')), 'V');

    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('N')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::dna15{}.assign_char('!')), 'N');
}

TEST(dna15, char_literal)
{
    EXPECT_EQ(seqan3::to_char('A'_dna15), 'A');
    EXPECT_EQ(seqan3::to_char('C'_dna15), 'C');
    EXPECT_EQ(seqan3::to_char('G'_dna15), 'G');

    EXPECT_EQ(seqan3::to_char('U'_dna15), 'T');
    EXPECT_EQ(seqan3::to_char('T'_dna15), 'T');

    EXPECT_EQ(seqan3::to_char('R'_dna15), 'R');
    EXPECT_EQ(seqan3::to_char('Y'_dna15), 'Y');
    EXPECT_EQ(seqan3::to_char('S'_dna15), 'S');
    EXPECT_EQ(seqan3::to_char('W'_dna15), 'W');
    EXPECT_EQ(seqan3::to_char('K'_dna15), 'K');
    EXPECT_EQ(seqan3::to_char('M'_dna15), 'M');
    EXPECT_EQ(seqan3::to_char('B'_dna15), 'B');
    EXPECT_EQ(seqan3::to_char('D'_dna15), 'D');
    EXPECT_EQ(seqan3::to_char('H'_dna15), 'H');
    EXPECT_EQ(seqan3::to_char('V'_dna15), 'V');

    EXPECT_EQ(seqan3::to_char('N'_dna15), 'N');
    EXPECT_EQ(seqan3::to_char('!'_dna15), 'N');
}

TEST(dna15, string_literal)
{
    seqan3::dna15_vector v;
    v.resize(5, 'A'_dna15);
    EXPECT_EQ(v, "AAAAA"_dna15);

    std::vector<seqan3::dna15> w{'A'_dna15, 'C'_dna15, 'G'_dna15, 'T'_dna15, 'U'_dna15, 'N'_dna15};
    EXPECT_EQ(w, "ACGTTN"_dna15);
}

TEST(dna15, char_is_valid)
{
    constexpr auto validator =
        seqan3::is_char<'A'> || seqan3::is_char<'C'> || seqan3::is_char<'G'> || seqan3::is_char<'T'>
        || seqan3::is_char<'U'> || seqan3::is_char<'a'> || seqan3::is_char<'c'> || seqan3::is_char<'g'>
        || seqan3::is_char<'t'> || seqan3::is_char<'u'> || seqan3::is_char<'N'> || seqan3::is_char<'n'>
        || seqan3::is_char<'R'> || seqan3::is_char<'Y'> || seqan3::is_char<'S'> || seqan3::is_char<'W'>
        || seqan3::is_char<'K'> || seqan3::is_char<'M'> || seqan3::is_char<'B'> || seqan3::is_char<'D'>
        || seqan3::is_char<'H'> || seqan3::is_char<'V'> || seqan3::is_char<'r'> || seqan3::is_char<'y'>
        || seqan3::is_char<'s'> || seqan3::is_char<'w'> || seqan3::is_char<'k'> || seqan3::is_char<'m'>
        || seqan3::is_char<'b'> || seqan3::is_char<'d'> || seqan3::is_char<'h'> || seqan3::is_char<'v'>;
    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(seqan3::dna15::char_is_valid(c), validator(c));
}
