// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <deque>

#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

namespace seqan3::detail
{
template <std::ranges::view urng_t>
class syncmer_view : public std::ranges::view_interface<syncmer_view<urng_t>>
{
private:
    static constexpr bool const_iterable = seqan3::const_iterable_range<urng_t>;

    urng_t urange{};
    size_t kmer_size{};
    size_t smer_size{};
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

    explicit syncmer_view(urng_t urange, size_t const kmer_size, size_t const smer_size, size_t const offset) :
        urange{std::move(urange)},
        kmer_size{kmer_size},
        smer_size{smer_size},
        offset{offset}
    {}

    template <typename other_urng_t>
        requires (std::ranges::viewable_range<other_urng_t>
                  && std::constructible_from<urng_t, std::ranges::ref_view<std::remove_reference_t<other_urng_t>>>)
    explicit syncmer_view(other_urng_t && urange, size_t const kmer_size, size_t const smer_size, size_t const offset) :
        urange{std::views::all(std::forward<other_urng_t>(urange))},
        kmer_size{kmer_size},
        smer_size{smer_size},
        offset{offset}
    {}

    template <typename other_urng_t>
        requires (std::ranges::viewable_range<other_urng_t>
                  && std::constructible_from<urng_t, std::views::all_t<other_urng_t>>)
    explicit syncmer_view(other_urng_t && urange, size_t const kmer_size, size_t const smer_size, size_t const offset) :
        urange{std::views::all(std::forward<other_urng_t>(urange))},
        kmer_size{kmer_size},
        smer_size{smer_size},
        offset{offset}
    {}

    basic_iterator<false> begin()
    {
        return {std::ranges::begin(urange), std::ranges::end(urange), kmer_size, smer_size, offset};
    }

    basic_iterator<true> begin() const
        requires const_iterable
    {
        return {std::ranges::begin(urange), std::ranges::end(urange), kmer_size, smer_size, offset};
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
    uint64_t smer_mask{};
    size_t offset{};

    value_type fwd_smer_value{};
    value_type rc_smer_value{};
    value_type fwd_kmer_value{};
    value_type rc_kmer_value{};
    value_type syncmer_value{};
    std::deque<value_type> fwd_smer_values{};
    std::deque<value_type> rc_smer_values{};
    size_t fwd_smer_position{};
    size_t rc_smer_position{};
    size_t kmer_size{};
    size_t smer_size{};

public:
    basic_iterator() = default;
    basic_iterator(basic_iterator const &) = default;
    basic_iterator(basic_iterator &&) = default;
    basic_iterator & operator=(basic_iterator const &) = default;
    basic_iterator & operator=(basic_iterator &&) = default;
    ~basic_iterator() = default;

    // basic_iterator(basic_iterator<!const_range> const & it)
    //     requires const_range
    //     :
    //     urng_iterator{std::move(it.urng_iterator)},
    //     urng_sentinel{std::move(it.urng_sentinel)},
    //     mask{std::move(it.mask)},
    //     offset{std::move(it.offset)},
    //     fwd_smer_value{std::move(it.fwd_smer_value)},
    //     fwd_kmer_value{std::move(it.fwd_kmer_value)},
    //     fwd_smer_values{std::move(it.fwd_smer_values)},
    //     fwd_smer_position{std::move(it.fwd_smer_position)}
    // {}

    basic_iterator(urng_iterator_t urng_iterator,
                   urng_sentinel_t urng_sentinel,
                   size_t const kmer_size,
                   size_t const smer_size,
                   size_t const offset) :
        urng_iterator{std::move(urng_iterator)},
        urng_sentinel{std::move(urng_sentinel)},
        mask{(1ULL << (2 * kmer_size)) - 1u},
        smer_mask{(1ULL << (2 * smer_size)) - 1u},
        offset{offset},
        kmer_size{kmer_size},
        smer_size{smer_size}
    {
        size_t const window_size = kmer_size - smer_size + 1u;

        init(window_size);
    }

    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return (lhs.urng_iterator == rhs.urng_iterator) && (rhs.urng2_iterator == rhs.urng2_iterator)
            && (lhs.fwd_smer_values.size() == rhs.fwd_smer_values.size());
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
        return syncmer_value;
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
        return offset;
    }

private:
    void update_kmer_value()
    {
        fwd_kmer_value <<= 2;
        fwd_kmer_value |= get_smer_value();
        fwd_kmer_value &= mask;

        rc_kmer_value >>= 2;
        rc_kmer_value |= get_smer_rc_value() << (2 * (kmer_size - smer_size));
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

    auto get_smer_rc_value() const
    {
        value_type current_smer = get_smer_value();
        value_type rc_value{};

        for (size_t i = 0u; i < smer_size - 1u; ++i)
        {
            rc_value |= current_smer & 3u;
            rc_value <<= 2;
            current_smer >>= 2;
        }

        rc_value |= current_smer & 3u;

        return rc_value ^ smer_mask;
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
            fwd_smer_values.push_back(get_smer_value());
            rc_smer_values.push_front(get_smer_rc_value());
            next_smer();
        }
        fwd_smer_values.push_back(get_smer_value());
        rc_smer_values.push_front(get_smer_rc_value());
        update_kmer_value();

        auto fwd_smer_it = std::ranges::min_element(fwd_smer_values, std::less_equal<value_type>{});
        fwd_smer_value = *fwd_smer_it;
        fwd_smer_position = std::distance(std::begin(fwd_smer_values), fwd_smer_it);

        auto rc_smer_it = std::ranges::min_element(rc_smer_values, std::less_equal<value_type>{});
        rc_smer_value = *rc_smer_it;
        rc_smer_position = std::distance(std::begin(rc_smer_values), rc_smer_it);

        if (fwd_kmer_value <= rc_kmer_value)
        {
            if (offset != fwd_smer_position)
                next_unique_syncmer();
            else
                syncmer_value = fwd_kmer_value;
        }
        else if (offset != rc_smer_position)
            next_unique_syncmer();
        else
            syncmer_value = rc_kmer_value;
    }

    bool next_syncmer()
    {
        next_smer();
        if (urng_iterator == urng_sentinel)
            return true;

        value_type const new_value = get_smer_value();
        value_type const new_rc_value = get_smer_rc_value();

        fwd_smer_values.pop_front();
        fwd_smer_values.push_back(new_value);
        update_kmer_value();

        rc_smer_values.pop_back();
        rc_smer_values.push_front(new_rc_value);

        if (fwd_smer_position == 0)
        {
            auto smer_it = std::ranges::min_element(fwd_smer_values, std::less_equal<value_type>{});
            fwd_smer_value = *smer_it;
            fwd_smer_position = std::distance(std::begin(fwd_smer_values), smer_it);
        }
        else if (new_value < fwd_smer_value)
        {
            fwd_smer_value = new_value;
            fwd_smer_position = fwd_smer_values.size() - 1;
        }
        else
        {
            --fwd_smer_position;
        }

        if (rc_smer_position == 0)
        {
            auto smer_it = std::ranges::min_element(rc_smer_values, std::less_equal<value_type>{});
            rc_smer_value = *smer_it;
            rc_smer_position = std::distance(std::begin(rc_smer_values), smer_it);
        }
        else if (new_rc_value < rc_smer_value)
        {
            rc_smer_value = new_rc_value;
            rc_smer_position = 0;
        }
        else
        {
            ++rc_smer_position;
        }

        if (fwd_kmer_value <= rc_kmer_value)
        {
            if (offset == fwd_smer_position)
            {
                syncmer_value = fwd_kmer_value;
                return true;
            }
        }
        else if (offset == rc_smer_position)
        {
            syncmer_value = rc_kmer_value;
            return true;
        }

        return false;
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
        return syncmer_view{view, kmer_size, smer_size, offset};
    }
};

} // namespace seqan3::detail

namespace seqan3::views
{

inline constexpr auto syncmer = detail::syncmer_fn{};

}
