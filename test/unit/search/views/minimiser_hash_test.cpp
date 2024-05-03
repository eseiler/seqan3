// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

// #include <gtest/gtest.h>

#include <bitset>
#include <forward_list>
#include <list>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_and_window_hash.hpp>
#include <seqan3/search/views/minimiser_hash_and_positions.hpp>
#include <seqan3/test/expect_range_eq.hpp>

// #include "../../range/iterator_test_template.hpp"

using seqan3::operator""_dna4;
// using seqan3::operator""_shape;
// using result_t = std::vector<size_t>;
// using iterator_type = std::ranges::iterator_t<
//     decltype(std::declval<seqan3::dna4_vector &>()
//              | seqan3::views::minimiser_hash(seqan3::ungapped{4}, seqan3::window_size{8}, seqan3::seed{0}))>;

// static constexpr seqan3::shape ungapped_shape = seqan3::ungapped{4};
// static constexpr seqan3::shape gapped_shape = 0b1001_shape;
// static constexpr auto ungapped_view =
//     seqan3::views::minimiser_hash(ungapped_shape, seqan3::window_size{8}, seqan3::seed{0});
// static constexpr auto gapped_view =
//     seqan3::views::minimiser_hash(gapped_shape, seqan3::window_size{8}, seqan3::seed{0});

// template <>
// struct iterator_fixture<iterator_type> : public ::testing::Test
// {
//     using iterator_tag = std::forward_iterator_tag;
//     static constexpr bool const_iterable = false;

//     seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
//     result_t expected_range{26, 97, 27, 6, 1};

//     using test_range_t = decltype(text | ungapped_view);
//     test_range_t test_range = text | ungapped_view;
// };

// using test_type = ::testing::Types<iterator_type>;
// INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_type, );

// template <typename T>
// class minimiser_hash_properties_test : public ::testing::Test
// {};

// using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
//                                                 std::vector<seqan3::dna4> const,
//                                                 seqan3::bitpacked_sequence<seqan3::dna4>,
//                                                 seqan3::bitpacked_sequence<seqan3::dna4> const,
//                                                 std::list<seqan3::dna4>,
//                                                 std::list<seqan3::dna4> const>;

// TYPED_TEST_SUITE(minimiser_hash_properties_test, underlying_range_types, );
// class minimiser_hash_test : public ::testing::Test
// {
// protected:
//     std::vector<seqan3::dna4> text1{"AAAAAAAAAAAAAAAAAAA"_dna4};
//     std::vector<seqan3::dna4> text1_short{"AAAAAA"_dna4};
//     result_t result1{0, 0, 0}; // Same for ungapped and gapped
//     result_t ungapped_default_seed{0x8F3F73B5CF1C9A21, 0x8F3F73B5CF1C9A21, 0x8F3F73B5CF1C9A21};
//     result_t gapped_default_seed{0x8F3F73B5CF1C9AD1, 0x8F3F73B5CF1C9AD1, 0x8F3F73B5CF1C9AD1};

//     std::vector<seqan3::dna4> text2{"AC"_dna4};
//     result_t result2{};

//     std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
//     result_t ungapped3{26, 97, 27, 6, 1}; // ACGG, CGAC, ACGT, aacg, aaac
//     result_t ungapped_stop_at_t3{26, 97}; // ACGG, CGAC
//     result_t gapped3{2, 5, 3, 2, 1};      // A--G, C--C, A--T, a--g, a--c "-" for gap
//     result_t gapped_stop_at_t3{2, 5};     // A--G, C--C "-" for gap
// };

// TEST(minimiser_pos, general)
// {

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

void minimiser_pos()
{
    // seqan3::dna4_vector text{"ACGTCGACGTTTAGAAAAAAAAAAAAAAAA"_dna4};
    // seqan3::dna4_vector text{"ACGTCGAC"_dna4};
    seqan3::dna4_vector text{"TCATCAGTAGCTACAATACG"_dna4};

    constexpr size_t const minimiser_size = 2;

    auto ungapped_view = bsc::views::minimiser_hash_and_positions({.minimiser_size = minimiser_size, .window_size = 5});

    // std::vector<size_t> ungapped_result{26, 97, 27, 6, 1}; // ACGG, CGAC, ACGT, aacg, aaac

    for (auto && res : text | ungapped_view)
    {
        seqan3::debug_stream << kmer_to_string(res.minimiser_value, minimiser_size) << ',' << res.range_position << ','
                             << res.occurrences << '\n';
    }

    // 27,0,0
    // 97,1,0
    // 27,2,4
    // 111,7,0
    // 32,8,1
    // 0,10,12

    // AT,0,2
    // AG,2,7
    // AC,9,2
    // AA,11,4
    // AC,15,1
}

void minimiser_both()
{
    // seqan3::dna4_vector text{"ACGTCGACGTTTAGAAAAAAAAAAAAAAAA"_dna4};
    // seqan3::dna4_vector text{"ACGTCGAC"_dna4};
    seqan3::dna4_vector text{"TCATCAGTAGCTACAATACG"_dna4};
    // seqan3::dna4_vector text{"TCAT"_dna4};

    constexpr size_t const minimiser_size = 2;
    constexpr size_t const window_size = 5;

    auto ungapped_view =
        bsc::views::minimiser_and_window_hash({.minimiser_size = minimiser_size, .window_size = window_size});

    // std::vector<size_t> ungapped_result{26, 97, 27, 6, 1}; // ACGG, CGAC, ACGT, aacg, aaac

    for (auto && res : text | ungapped_view)
    {
        seqan3::debug_stream << kmer_to_string(res.minimiser_value, minimiser_size) << ','
                             << kmer_to_string(res.window_value, window_size) << '\n';
    }
}

int main()
{
    minimiser_pos();
    // minimiser_both();
}
// }

// 0  ACGTCGAC | GTTTAGAAAAAAAAAAAAAAAA
// 1  CGTCGACG | TTTAGAAAAAAAAAAAAAAAA
// 2  GTCGACGT | TTAGAAAAAAAAAAAAAAAA
// 3  TCGACGTT | TAGAAAAAAAAAAAAAAAA
// 4  CGACGTTT | AGAAAAAAAAAAAAAAAA
// 5  GACGTTTA | GAAAAAAAAAAAAAAAA
// 6  ACGTTTAG | AAAAAAAAAAAAAAAA
// 7  CGTTTAGA | AAAAAAAAAAAAAAA
// 8  GTTTAGAA | AAAAAAAAAAAAAA
// 9  TTTAGAAA | AAAAAAAAAAAAA
// 0  TTAGAAAA | AAAAAAAAAAAA
// 1  TAGAAAAA | AAAAAAAAAAA
// 2  AGAAAAAA | AAAAAAAAAA
// 3  GAAAAAAA | AAAAAAAAA
// 4  AAAAAAAA | AAAAAAAA
// 5  AAAAAAAA | AAAAAAA
// 6  AAAAAAAA | AAAAAA
// 7  AAAAAAAA | AAAAA
// 8  AAAAAAAA | AAAA
// 9  AAAAAAAA | AAA
// 0  AAAAAAAA | AA
// 1  AAAAAAAA | A
// 2  AAAAAAAA |

// TYPED_TEST(minimiser_hash_properties_test, different_input_ranges)
// {
//     TypeParam text{'A'_dna4,
//                    'C'_dna4,
//                    'G'_dna4,
//                    'T'_dna4,
//                    'C'_dna4,
//                    'G'_dna4,
//                    'A'_dna4,
//                    'C'_dna4,
//                    'G'_dna4,
//                    'T'_dna4,
//                    'T'_dna4,
//                    'T'_dna4,
//                    'A'_dna4,
//                    'G'_dna4};            // ACGTCGACGTTTAG
//     result_t ungapped{27, 97, 27, 6, 1}; // ACGT, CGAC, ACGT, aacg, aaac
//     result_t gapped{3, 5, 3, 2, 1};      // A--T, C--C, A--T, a--g, a--c - "-" for gap
//     EXPECT_RANGE_EQ(ungapped, text | ungapped_view);
//     EXPECT_RANGE_EQ(gapped, text | gapped_view);
// }

// TEST_F(minimiser_hash_test, ungapped)
// {
//     EXPECT_RANGE_EQ(result1, text1 | ungapped_view);
//     EXPECT_RANGE_EQ(result2, text2 | ungapped_view);
//     EXPECT_RANGE_EQ(ungapped3, text3 | ungapped_view);

//     auto stop_at_t = std::views::take_while(
//         [](seqan3::dna4 const x)
//         {
//             return x != 'T'_dna4;
//         });
//     EXPECT_RANGE_EQ(ungapped_stop_at_t3, text3 | stop_at_t | ungapped_view);
// }

// TEST_F(minimiser_hash_test, gapped)
// {
//     EXPECT_RANGE_EQ(result1, text1 | gapped_view);
//     EXPECT_RANGE_EQ(result2, text2 | gapped_view);
//     EXPECT_RANGE_EQ(gapped3, text3 | gapped_view);

//     auto stop_at_t = std::views::take_while(
//         [](seqan3::dna4 const x)
//         {
//             return x != 'T'_dna4;
//         });
//     EXPECT_RANGE_EQ(gapped_stop_at_t3, text3 | stop_at_t | gapped_view);
// }

// TEST_F(minimiser_hash_test, seed)
// {
//     EXPECT_RANGE_EQ(ungapped_default_seed,
//                     text1 | seqan3::views::minimiser_hash(ungapped_shape, seqan3::window_size{8}));
//     EXPECT_RANGE_EQ(gapped_default_seed, text1 | seqan3::views::minimiser_hash(gapped_shape, seqan3::window_size{8}));
// }

// TEST_F(minimiser_hash_test, shape_bigger_than_window)
// {
//     EXPECT_THROW(text1 | seqan3::views::minimiser_hash(ungapped_shape, seqan3::window_size{3}), std::invalid_argument);
//     EXPECT_THROW(text1 | seqan3::views::minimiser_hash(gapped_shape, seqan3::window_size{3}), std::invalid_argument);
// }
