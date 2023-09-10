// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <deque>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/utility/range/concept.hpp>

namespace seqan3::detail
{

struct syncmer_params
{
    size_t kmer_size{};
    size_t smer_size{};
    size_t offset{};
};

template <std::ranges::view urng_t>
class syncmer_view : public std::ranges::view_interface<syncmer_view<urng_t>>
{
private:
    urng_t urange{};
    syncmer_params params{};

    template <bool const_range>
    class basic_iterator;

public:
    syncmer_view()
        requires std::default_initializable<urng_t>
    = default;
    syncmer_view(syncmer_view const & rhs) = default;
    syncmer_view(syncmer_view && rhs) = default;
    syncmer_view & operator=(syncmer_view const & rhs) = default;
    syncmer_view & operator=(syncmer_view && rhs) = default;
    ~syncmer_view() = default;

    explicit syncmer_view(urng_t urange, syncmer_params const & params) : urange{std::move(urange)}, params{params}
    {}

    template <std::ranges::viewable_range other_urng_t>
        requires std::constructible_from<urng_t, std::views::all_t<other_urng_t>>
    explicit syncmer_view(other_urng_t && urange, syncmer_params const & params) :
        urange{std::views::all(std::forward<other_urng_t>(urange))},
        params{params}
    {}

    basic_iterator<false> begin()
    {
        return {std::ranges::begin(urange), std::ranges::end(urange), params};
    }

    basic_iterator<true> begin() const
        requires const_iterable_range<urng_t>
    {
        return {std::ranges::begin(urange), std::ranges::end(urange), params};
    }

    auto end() noexcept
    {
        return std::ranges::end(urange);
    }

    auto end() const noexcept
        requires const_iterable_range<urng_t>
    {
        return std::ranges::cend(urange);
    }
};

template <std::ranges::view urng_t>
template <bool const_range>
class syncmer_view<urng_t>::basic_iterator
{
public:
    using difference_type = std::ranges::range_difference_t<urng_t>;
    using value_type = uint64_t;
    using pointer = void;
    using reference = value_type;
    using iterator_category = std::forward_iterator_tag;
    using iterator_concept = iterator_category;

private:
    template <bool>
    friend class basic_iterator;

    using urng_iterator_t = maybe_const_iterator_t<const_range, urng_t>;
    using urng_sentinel_t = maybe_const_sentinel_t<const_range, urng_t>;

    urng_iterator_t text_it{};
    urng_sentinel_t text_end{};
    syncmer_params params{};
    size_t fwd_kmer_mask{};
    size_t fwd_smer_mask{};
    size_t rc_kmer_shift{};
    size_t rc_smer_shift{};

    value_type fwd_min_smer_value{};
    value_type rc_min_smer_value{};
    value_type fwd_kmer_value{};
    value_type rc_kmer_value{};
    value_type syncmer_value{};
    size_t fwd_smer_position{};
    size_t rc_smer_position{};
    std::deque<value_type> fwd_smer_values{};
    std::deque<value_type> rc_smer_values{};

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
        text_it{std::move(it.text_it)},
        text_end{std::move(it.text_end)},
        params{std::move(it.params)},
        fwd_kmer_mask{std::move(it.fwd_kmer_mask)},
        fwd_smer_mask{std::move(it.fwd_smer_mask)},
        rc_kmer_shift{std::move(it.rc_kmer_shift)},
        rc_smer_shift{std::move(it.rc_smer_shift)},
        fwd_min_smer_value{std::move(it.fwd_min_smer_value)},
        rc_min_smer_value{std::move(it.rc_min_smer_value)},
        fwd_kmer_value{std::move(it.fwd_kmer_value)},
        rc_kmer_value{std::move(it.rc_kmer_value)},
        fwd_smer_position{std::move(it.fwd_smer_position)},
        rc_smer_position{std::move(it.rc_smer_position)},
        fwd_smer_values{std::move(it.fwd_smer_values)},
        rc_smer_values{std::move(it.rc_smer_values)}
    {}

    basic_iterator(urng_iterator_t urng_iterator, urng_sentinel_t urng_sentinel, syncmer_params const & params) :
        text_it{std::move(urng_iterator)},
        text_end{std::move(urng_sentinel)},
        params{params},
        fwd_kmer_mask{(1ULL << (2 * params.kmer_size)) - 1u}, // k = 32?
        fwd_smer_mask{(1ULL << (2 * params.smer_size)) - 1u},
        rc_kmer_shift{2 * (params.kmer_size - 1u)},
        rc_smer_shift{2 * (params.smer_size - 1u)}
    {
        init();
    }

    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return (lhs.text_it == rhs.text_it);
    }

    friend bool operator!=(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return !(lhs == rhs);
    }

    friend bool operator==(basic_iterator const & lhs, urng_sentinel_t const & rhs)
    {
        return lhs.text_it == rhs;
    }

    friend bool operator==(urng_sentinel_t const & lhs, basic_iterator const & rhs)
    {
        return rhs == lhs;
    }

    friend bool operator!=(urng_sentinel_t const & lhs, basic_iterator const & rhs)
    {
        return !(lhs == rhs);
    }

    friend bool operator!=(basic_iterator const & lhs, urng_sentinel_t const & rhs)
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
        return text_it;
    }

    constexpr urng_iterator_t base() &&
    {
        return std::move(text_it);
    }

private:
    template <bool with_pop = true>
    void update_values()
    {
        value_type const new_rank = to_rank(*text_it);

        fwd_kmer_value <<= 2;
        fwd_kmer_value |= new_rank;
        fwd_kmer_value &= fwd_kmer_mask;

        value_type new_fwd_smer_value{fwd_smer_values.back()};
        new_fwd_smer_value <<= 2;
        new_fwd_smer_value |= new_rank;
        new_fwd_smer_value &= fwd_smer_mask;

        if constexpr (with_pop)
            fwd_smer_values.pop_front();
        fwd_smer_values.push_back(new_fwd_smer_value);

        rc_kmer_value >>= 2;
        rc_kmer_value |= (new_rank ^ 3u) << rc_kmer_shift;

        value_type new_rc_smer_value{rc_smer_values.front()};
        new_rc_smer_value >>= 2;
        new_rc_smer_value |= (new_rank ^ 3u) << rc_smer_shift;

        if constexpr (with_pop)
            rc_smer_values.pop_back();
        rc_smer_values.push_front(new_rc_smer_value);
    }

    void next_unique_syncmer()
    {
        while (!next_syncmer())
        {}
    }

    void find_minimum_fwd_smer()
    {
        auto fwd_smer_it = std::ranges::min_element(fwd_smer_values, std::less_equal<value_type>{});
        fwd_min_smer_value = *fwd_smer_it;
        fwd_smer_position = std::distance(std::begin(fwd_smer_values), fwd_smer_it);
    }

    void find_minimum_rc_smer()
    {
        auto rc_smer_it = std::ranges::min_element(rc_smer_values, std::less_equal<value_type>{});
        rc_min_smer_value = *rc_smer_it;
        rc_smer_position = std::distance(std::begin(rc_smer_values), rc_smer_it);
    }

    void init()
    {
        // Initial values for update_values()
        fwd_smer_values.push_back(value_type{});
        rc_smer_values.push_front(value_type{});

        // Don't keep smer values in queue while processing the first s-1 chracters.
        for (size_t i = 0u; i < params.smer_size - 1; ++i)
        {
            update_values();
            ++text_it;
        }
        // Fill queue with smer values.
        for (size_t i = params.smer_size - 1; i < params.kmer_size - 1u; ++i)
        {
            update_values<false>();
            ++text_it;
        }
        update_values<false>();

        // Remove initial values for update_values()
        fwd_smer_values.pop_front();
        rc_smer_values.pop_back();

        find_minimum_fwd_smer();
        find_minimum_rc_smer();

        if (fwd_kmer_value <= rc_kmer_value)
        {
            if (params.offset != fwd_smer_position)
                next_unique_syncmer();
            else
                syncmer_value = fwd_kmer_value;
        }
        else if (params.offset != rc_smer_position)
            next_unique_syncmer();
        else
            syncmer_value = rc_kmer_value;
    }

    bool next_syncmer()
    {
        ++text_it;

        if (text_it == text_end)
            return true;

        update_values();

        if (fwd_smer_position == 0)
        {
            find_minimum_fwd_smer();
        }
        else if (fwd_smer_values.back() < fwd_min_smer_value)
        {
            fwd_min_smer_value = fwd_smer_values.back();
            fwd_smer_position = fwd_smer_values.size() - 1u;
        }
        else
        {
            --fwd_smer_position;
        }

        if (rc_smer_position == 0)
        {
            find_minimum_rc_smer();
        }
        else if (rc_smer_values.front() < rc_min_smer_value)
        {
            rc_min_smer_value = rc_smer_values.front();
            rc_smer_position = 0;
        }
        else
        {
            ++rc_smer_position;
        }

        if (fwd_kmer_value <= rc_kmer_value)
        {
            if (params.offset == fwd_smer_position)
            {
                syncmer_value = fwd_kmer_value;
                return true;
            }
        }
        else if (params.offset == rc_smer_position)
        {
            syncmer_value = rc_kmer_value;
            return true;
        }

        return false;
    }
};

template <std::ranges::viewable_range rng_t>
syncmer_view(rng_t &&, syncmer_params const & params) -> syncmer_view<std::views::all_t<rng_t>>;

struct syncmer_fn
{
    constexpr auto operator()(syncmer_params const & params) const
    {
        return adaptor_from_functor{*this, params};
    }

    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, syncmer_params const & params) const
    {
        static_assert(std::same_as<std::ranges::range_value_t<urng_t>, seqan3::dna4>, "Only dna4 supported.");
        if (params.kmer_size == 0u)
            throw std::invalid_argument{"kmer_size must be > 0."};
        if (params.smer_size == 0u)
            throw std::invalid_argument{"smer_size must be > 0."};
        if (params.kmer_size < params.smer_size)
            throw std::invalid_argument{"kmer_size must be >= smer_size."};
        if (params.offset > params.kmer_size - params.smer_size)
            throw std::invalid_argument{"offset must be in [0, kmer_size - smer_size]."};
        return syncmer_view{std::forward<urng_t>(urange), params};
    }
};

} // namespace seqan3::detail

namespace seqan3::views
{

inline constexpr auto syncmer = detail::syncmer_fn{};

}
