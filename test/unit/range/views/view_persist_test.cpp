// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/unique.hpp>

#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

using namespace seqan3;

// ============================================================================
//  test templates
// ============================================================================

TEST(view_persist, delegate_to_view_all)
{
    std::string vec{"foo"};

    // pipe notation
    auto v = vec | views::persist;
    EXPECT_EQ("foo", std::string{v});

    // function notation
    std::string v2 = views::persist(vec) | std::ranges::to<std::string>;
    EXPECT_EQ("foo", v2);

    // combinability
    auto v3 = vec | views::persist | ranges::view::unique;
    EXPECT_EQ("fo", std::string{v3});
    std::string v3b = vec | std::views::reverse | views::persist | ranges::view::unique | std::ranges::to<std::string>;
    EXPECT_EQ("of", v3b);

    // store combined
    auto a1 = views::persist | ranges::view::unique;
    auto v5 = vec | a1;
    EXPECT_EQ("fo", std::string{v5});
}

TEST(view_persist, wrap_temporary)
{
    // pipe notation
    auto v = std::string{"foo"} | views::persist;
    EXPECT_EQ("foo", std::string(v));

    // function notation
    std::string v2 = views::persist(std::string{"foo"}) | std::ranges::to<std::string>;
    EXPECT_EQ("foo", v2);

    // combinability
    auto v3 = std::string{"foo"} | views::persist | ranges::view::unique;
    EXPECT_EQ("fo", std::string(v3));
    std::string v3b = std::string{"foo"}
                    | views::persist
                    | std::views::filter(is_char<'o'>)
                    | ranges::view::unique
                    | std::ranges::to<std::string>;
    EXPECT_EQ("o", v3b);
}

TEST(view_persist, const_)
{
    // inner const
    using t = std::string const;
    auto v = t{"foo"} | views::persist;
    EXPECT_EQ("foo", std::string{v});

    // outer const
    auto const & v2 = std::string{"foo"} | views::persist;
    EXPECT_EQ("foo", std::string{v2});

    // inner + outer const
    using t = std::string const;
    auto const & v3 = t{"foo"} | views::persist;
    EXPECT_EQ("foo", std::string{v3});
}

TEST(view_persist, concepts)
{
    std::string vec{"foobar"};
    EXPECT_TRUE(std::ranges::input_range<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(std::string{"foo"})>);
    EXPECT_FALSE(std::ranges::view<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::common_range<decltype(std::string{"foo"})>);
    EXPECT_TRUE(const_iterable_range<decltype(std::string{"foo"})>);
    EXPECT_TRUE((std::ranges::output_range<decltype(std::string{"foo"}), char>));

    auto v1 = std::string{"foo"} | views::persist;

    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(const_iterable_range<decltype(v1)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), char>));
}