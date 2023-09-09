// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <deque>

#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/utility/range/concept.hpp>
#include <seqan3/utility/type_traits/lazy_conditional.hpp>

namespace seqan3::detail
{
template <std::ranges::view urng_t>
class syncmer_view : public std::ranges::view_interface<syncmer_view<urng_t>>
{
private:
    static constexpr bool const_iterable = seqan3::const_iterable_range<urng_t>;

    urng_t urange{};
    size_t window_size{};
    size_t offset{};

    template <bool const_range>
    class basic_iterator;

    using sentinel = std::default_sentinel_t;

public:
    syncmer_view()
        requires std::default_initializable<urng_t>
    = default;
    syncmer_view(syncmer_view const & rhs) = default;
    syncmer_view(syncmer_view && rhs) = default;
    syncmer_view & operator=(syncmer_view const & rhs) = default;
    syncmer_view & operator=(syncmer_view && rhs) = default;
    ~syncmer_view() = default;

    explicit syncmer_view(urng_t urange, size_t const window_size, size_t const offset) :
        urange{std::move(urange)},
        window_size{window_size},
        offset{offset}
    {}

    template <typename other_urng_t>
        requires (std::ranges::viewable_range<other_urng_t>
                  && std::constructible_from<urng_t, std::ranges::ref_view<std::remove_reference_t<other_urng_t>>>)
    explicit syncmer_view(other_urng_t && urange, size_t const window_size, size_t const offset) :
        urange{std::views::all(std::forward<other_urng_t>(urange))},
        window_size{window_size},
        offset{offset}
    {}

    template <typename other_urng_t>
        requires (std::ranges::viewable_range<other_urng_t>
                  && std::constructible_from<urng_t, std::views::all_t<other_urng_t>>)
    explicit syncmer_view(other_urng_t && urange, size_t const window_size, size_t const offset) :
        urange{std::views::all(std::forward<other_urng_t>(urange))},
        window_size{window_size},
        offset{offset}
    {}

    basic_iterator<false> begin()
    {
        return {std::ranges::begin(urange), std::ranges::end(urange), window_size, offset};
    }

    basic_iterator<true> begin() const
        requires const_iterable
    {
        return {std::ranges::begin(urange), std::ranges::end(urange), window_size, offset};
    }

    sentinel end() const
    {
        return {};
    }
};

template <std::ranges::view urng_t>
template <bool const_range>
class syncmer_view<urng_t>::basic_iterator
{
public:
    using difference_type = std::ranges::range_difference_t<urng_t>;
    using value_type = std::ranges::range_value_t<urng_t>;
    using pointer = void;
    using reference = value_type;
    using iterator_category = std::forward_iterator_tag;
    using iterator_concept = iterator_category;

private:
    template <bool>
    friend class basic_iterator;

    using urng_sentinel_t = maybe_const_sentinel_t<const_range, urng_t>;
    using urng_iterator_t = maybe_const_iterator_t<const_range, urng_t>;

    urng_iterator_t urng_iterator{};
    urng_sentinel_t urng_sentinel{};
    uint64_t mask{};
    size_t offset{};

    value_type smer_value{};
    value_type kmer_value{};
    std::deque<value_type> smer_values{};
    size_t smer_position_offset{};

public:
    basic_iterator() = default;
    basic_iterator(basic_iterator const &) = default;
    basic_iterator(basic_iterator &&) = default;
    basic_iterator & operator=(basic_iterator const &) = default;
    basic_iterator & operator=(basic_iterator &&) = default;
    ~basic_iterator() = default;

    basic_iterator(basic_iterator<!const_range> const & it)
        requires const_range
        :
        urng_iterator{std::move(it.urng_iterator)},
        urng_sentinel{std::move(it.urng_sentinel)},
        mask{std::move(it.mask)},
        offset{std::move(it.offset)},
        smer_value{std::move(it.smer_value)},
        kmer_value{std::move(it.kmer_value)},
        smer_values{std::move(it.smer_values)},
        smer_position_offset{std::move(it.smer_position_offset)}
    {}

    basic_iterator(urng_iterator_t urng_iterator,
                   urng_sentinel_t urng_sentinel,
                   size_t window_size,
                   size_t const offset) :
        urng_iterator{std::move(urng_iterator)},
        urng_sentinel{std::move(urng_sentinel)},
        offset{offset}
    {
        size_t size = std::ranges::distance(urng_iterator, urng_sentinel);
        window_size = std::min<size_t>(window_size, size);
        mask = (1ULL << (2 * window_size)) - 1u;

        init(window_size);
    }

    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return (lhs.urng_iterator == rhs.urng_iterator) && (rhs.urng2_iterator == rhs.urng2_iterator)
            && (lhs.smer_values.size() == rhs.smer_values.size());
    }

    friend bool operator!=(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return !(lhs == rhs);
    }

    friend bool operator==(basic_iterator const & lhs, sentinel const &)
    {
        return lhs.urng_iterator == lhs.urng_sentinel;
    }

    friend bool operator==(sentinel const & lhs, basic_iterator const & rhs)
    {
        return rhs == lhs;
    }

    friend bool operator!=(sentinel const & lhs, basic_iterator const & rhs)
    {
        return !(lhs == rhs);
    }

    friend bool operator!=(basic_iterator const & lhs, sentinel const & rhs)
    {
        return !(lhs == rhs);
    }

    basic_iterator & operator++() noexcept
    {
        next_unique_syncmer();
        return *this;
    }

    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        next_unique_syncmer();
        return tmp;
    }

    value_type operator*() const noexcept
    {
        return kmer_value;
    }

    constexpr urng_iterator_t const & base() const & noexcept
    {
        return urng_iterator;
    }

    constexpr urng_iterator_t base() &&
    {
        return std::move(urng_iterator);
    }

    size_t get_offset() const noexcept
    {
        return smer_position_offset;
    }

private:
    void update_kmer_value()
    {
        kmer_value <<= 2;
        kmer_value |= get_smer_value();
        kmer_value &= mask;
    }

    void next_unique_syncmer()
    {
        while (!next_syncmer())
        {}
    }

    auto get_smer_value() const
    {
        return *urng_iterator;
    }

    void next_smer()
    {
        ++urng_iterator;
    }

    void init(size_t const window_size)
    {
        if (window_size == 0u)
            return;

        for (size_t i = 0u; i < window_size - 1u; ++i)
        {
            update_kmer_value();
            smer_values.push_back(get_smer_value());
            next_smer();
        }
        smer_values.push_back(get_smer_value());
        update_kmer_value();
        auto smer_it = std::ranges::min_element(smer_values, std::less_equal<value_type>{});
        smer_value = *smer_it;
        smer_position_offset = std::distance(std::begin(smer_values), smer_it);
        if (offset != smer_position_offset)
            next_unique_syncmer();
    }

    bool next_syncmer()
    {
        next_smer();
        if (urng_iterator == urng_sentinel)
            return true;

        value_type const new_value = get_smer_value();

        smer_values.pop_front();
        smer_values.push_back(new_value);
        update_kmer_value();

        if (smer_position_offset == 0)
        {
            auto smer_it = std::ranges::min_element(smer_values, std::less_equal<value_type>{});
            smer_value = *smer_it;
            smer_position_offset = std::distance(std::begin(smer_values), smer_it);
            return smer_position_offset == offset;
        }

        if (new_value < smer_value)
        {
            smer_value = new_value;
            smer_position_offset = smer_values.size() - 1;
            return smer_position_offset == offset;
        }

        --smer_position_offset;
        return smer_position_offset == offset;
    }
};

template <std::ranges::viewable_range rng_t>
syncmer_view(rng_t &&, size_t const window_size, size_t const offset) -> syncmer_view<std::views::all_t<rng_t>>;

struct syncmer_fn
{
    constexpr auto operator()(size_t const kmer_size, uint8_t const smer_size, size_t const offset) const
    {
        return adaptor_from_functor{*this, kmer_size, smer_size, offset};
    }

    template <std::ranges::range urng_t>
    constexpr auto
    operator()(urng_t && urange, size_t const kmer_size, uint8_t const smer_size, size_t const offset) const
    {
        auto view = std::forward<urng_t>(urange) | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{smer_size}});
        size_t const window_size = kmer_size - smer_size + 1u;

        return syncmer_view{view, window_size, offset};
    }
};

} // namespace seqan3::detail

namespace seqan3::views
{

inline constexpr auto syncmer = detail::syncmer_fn{};

}
