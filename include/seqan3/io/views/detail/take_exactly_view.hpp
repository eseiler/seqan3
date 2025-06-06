// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::views::take_exactly and seqan3::views::take_exactly_or_throw.
 */

#pragma once

#include <algorithm>
#include <concepts>
#include <iterator>
#include <ranges>
#include <span>
#include <type_traits>

#include <seqan3/core/detail/iterator_traits.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/core/range/detail/inherited_iterator_base.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/utility/range/concept.hpp>
#include <seqan3/utility/type_traits/detail/transformation_trait_or.hpp>

namespace seqan3::detail
{

// ============================================================================
//  view_take_exactly
// ============================================================================

/*!\brief The type returned by seqan3::views::take, seqan3::detail::take_exactly and
 *        seqan3::detail::take_exactly_or_throw.
 * \tparam urng_t    The type of the underlying ranges, must satisfy seqan3::views::concept.
 * \tparam or_throw  Whether to throw an exception when the input is exhausted before the end of line is reached.
 * \implements std::ranges::view
 * \implements std::ranges::random_access_range
 * \implements std::ranges::sized_range
 * \ingroup io_views
 *
 * \details
 *
 * Note that most members of this class are generated by std::ranges::view_interface which is not yet documented here.
 */
template <std::ranges::view urng_t, bool or_throw>
class view_take_exactly : public std::ranges::view_interface<view_take_exactly<urng_t, or_throw>>
{
private:
    //!\brief The underlying range.
    urng_t urange;

    //!\brief The desired target_size.
    size_t target_size;

    //!\brief The forward declared iterator type.
    template <bool const_range>
    class basic_iterator;

private:
    /*!\name Associated types
     * \{
     */
    //!\brief The iterator type of this view (a random access iterator).
    using iterator = basic_iterator<false>;
    /*!\brief Note that this declaration does not give any compiler errors for non-const iterable ranges. Although
     * `basic_iterator` inherits from std::ranges::iterator_t which is not defined on a const-range, i.e. `urng_t const,
     *  if it is not const-iterable. We only just declare this type and never instantiate it, i.e. use this type within
     *  this class, if the underlying range is not const-iterable.
     */
    using const_iterator = basic_iterator<true>;
    //!\}

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    view_take_exactly() = default;                                          //!< Defaulted.
    view_take_exactly(view_take_exactly const & rhs) = default;             //!< Defaulted.
    view_take_exactly(view_take_exactly && rhs) = default;                  //!< Defaulted.
    view_take_exactly & operator=(view_take_exactly const & rhs) = default; //!< Defaulted.
    view_take_exactly & operator=(view_take_exactly && rhs) = default;      //!< Defaulted.
    ~view_take_exactly() = default;                                         //!< Defaulted.

    /*!\brief Construct from another View.
     * \param[in] _urange The underlying range.
     * \param[in] _size   The desired size (after which to stop returning elements).
     * \throws unexpected_end_of_input If `or_throw && seqan3::sized_range<urng_t>`.
     */
    constexpr view_take_exactly(urng_t _urange, size_t const _size) : urange{std::move(_urange)}, target_size{_size}
    {
        if constexpr (std::ranges::sized_range<urng_t>)
        {
            if (std::ranges::size(urange) < target_size)
            {
                if constexpr (or_throw)
                {
                    throw std::invalid_argument{
                        "You are trying to construct a detail::take_exactly_or_throw from a range that is strictly "
                        "smaller."};
                }
                else
                {
                    target_size = std::ranges::size(urange);
                }
            }
        }
    }

    /*!\brief Construct from another viewable_range.
     * \tparam rng_t      Type of the passed range; `urng_t` must be constructible from this.
     * \param[in] _urange The underlying range.
     * \param[in] _size   The desired size (after which to stop returning elements).
     * \throws unexpected_end_of_input If `or_throw && seqan3::sized_range<urng_t>`.
     */
    template <std::ranges::viewable_range rng_t>
        requires std::constructible_from<rng_t, std::views::all_t<rng_t>>
    constexpr view_take_exactly(rng_t && _urange, size_t const _size) :
        view_take_exactly{std::views::all(std::forward<rng_t>(_urange)), _size}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to end().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    constexpr auto begin() noexcept
    {
        if constexpr (std::ranges::random_access_range<urng_t> && std::ranges::sized_range<urng_t>)
            return std::ranges::begin(urange);
        else
            return iterator{std::ranges::begin(urange), 0, target_size, this};
    }

    //!\copydoc begin()
    constexpr auto begin() const noexcept
        requires const_iterable_range<urng_t> && std::ranges::forward_range<urng_t>
    {
        if constexpr (std::ranges::random_access_range<urng_t> && std::ranges::sized_range<urng_t>)
            return std::ranges::cbegin(urange);
        else
        {
            // const_iterator does not work if the underlying iterator is a std::input_iterator (it needs to access
            // target_size)
            return const_iterator{std::ranges::begin(urange), 0, target_size};
        }
    }

    /*!\brief Returns an iterator to the element following the last element of the range.
     * \returns Iterator to the end.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    constexpr auto end() noexcept
    {
        if constexpr (std::ranges::random_access_range<urng_t> && std::ranges::sized_range<urng_t>)
            return std::ranges::begin(urange) + target_size;
        else
            return std::ranges::end(urange);
    }

    //!\copydoc end()
    constexpr auto end() const noexcept
        requires const_iterable_range<urng_t> && std::ranges::forward_range<urng_t>
    {
        if constexpr (std::ranges::random_access_range<urng_t> && std::ranges::sized_range<urng_t>)
            return std::ranges::cbegin(urange) + target_size;
        else
            return std::ranges::cend(urange);
    }
    //!\}

    /*!\brief Returns the number of elements in the view.
     * \returns The number of elements in the view.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    constexpr auto size() const noexcept
    {
        return target_size;
    }
};

//!\brief Template argument type deduction guide that strips references.
//!\relates seqan3::detail::view_take_exactly
template <typename urng_t, bool or_throw = false>
view_take_exactly(urng_t &&, size_t) -> view_take_exactly<std::views::all_t<urng_t>, or_throw>;

//!\brief The iterator for the view_take_exactly. It inherits from the underlying type, but overwrites several operators.
//!\tparam rng_t Should be `urng_t` for defining #iterator and `urng_t const` for defining #const_iterator.
template <std::ranges::view urng_t, bool or_throw>
template <bool const_range>
class view_take_exactly<urng_t, or_throw>::basic_iterator :
    public inherited_iterator_base<basic_iterator<const_range>, maybe_const_iterator_t<const_range, urng_t>>
{
private:
    //!\brief The iterator type of the underlying range.
    using base_base_t = maybe_const_iterator_t<const_range, urng_t>;
    //!\brief The CRTP wrapper type.
    using base_t = inherited_iterator_base<basic_iterator, maybe_const_iterator_t<const_range, urng_t>>;

    //!\brief The sentinel type is identical to that of the underlying range.
    using sentinel_type = maybe_const_sentinel_t<const_range, urng_t>;

    //!\brief The current position.
    size_t pos{};

    //!\brief The size parameter to the view.
    size_t max_pos{};

    //!\brief A pointer to host, s.t. the size of the view can shrink on pure input ranges.
    std::conditional_t<!std::forward_iterator<base_base_t>, view_take_exactly *, detail::ignore_t> host_ptr;

public:
    /*!\name Constructors, destructor and assignment
     * \brief Exceptions specification is implicitly inherited.
     * \{
     */
    basic_iterator() = default;                                       //!< Defaulted.
    basic_iterator(basic_iterator const & rhs) = default;             //!< Defaulted.
    basic_iterator(basic_iterator && rhs) = default;                  //!< Defaulted.
    basic_iterator & operator=(basic_iterator const & rhs) = default; //!< Defaulted.
    basic_iterator & operator=(basic_iterator && rhs) = default;      //!< Defaulted.
    ~basic_iterator() = default;                                      //!< Defaulted.

    //!\brief Constructor that delegates to the CRTP layer.
    constexpr basic_iterator(base_base_t it) noexcept(noexcept(base_t{it})) : base_t{std::move(it)}
    {}

    //!\brief Constructor that delegates to the CRTP layer and initialises the members.
    constexpr basic_iterator(base_base_t it,
                             size_t const _pos,
                             size_t const _max_pos,
                             view_take_exactly * host = nullptr) noexcept(noexcept(base_t{it})) :
        base_t{std::move(it)},
        pos{_pos},
        max_pos(_max_pos)
    {
        host_ptr = host;

        if constexpr (!std::forward_iterator<base_base_t>)
        {
            assert(host_ptr != nullptr);
        }
    }
    //!\}

    using typename base_t::difference_type;
    using typename base_t::reference;

    /*!\name Arithmetic operators
     * \brief seqan3::detail::inherited_iterator_base operators are used unless specialised here.
     * \{
     */

    //!\brief Increments the iterator by one.
    constexpr basic_iterator & operator++() noexcept(noexcept(++std::declval<base_t &>()))
    {
        base_t::operator++();
        ++pos;
        if constexpr (!std::forward_iterator<base_base_t>)
            --host_ptr->target_size;
        return *this;
    }

    //!\brief Returns an iterator incremented by one.
    constexpr decltype(auto) operator++(int) noexcept(noexcept(++std::declval<basic_iterator &>())
                                                      && (std::same_as<decltype(std::declval<base_base_t &>()++), void>
                                                          || std::is_nothrow_copy_constructible_v<basic_iterator>))
    {
        // if underlying iterator is a C++20 input iterator (i.e. returns void), return void too.
        if constexpr (std::same_as<decltype(std::declval<base_base_t &>()++), void>)
        {
            ++(*this);
        }
        else
        {
            basic_iterator cpy{*this};
            ++(*this);
            return cpy;
        }
    }

    //!\brief Decrements the iterator by one.
    constexpr basic_iterator & operator--() noexcept(noexcept(--std::declval<base_base_t &>()))
        requires std::bidirectional_iterator<base_base_t>
    {
        base_t::operator--();
        --pos;
        return *this;
    }

    //!\brief Returns an iterator decremented by one.
    constexpr basic_iterator operator--(int)
        noexcept(noexcept(--std::declval<basic_iterator &>()) && std::is_nothrow_copy_constructible_v<basic_iterator>)
        requires std::bidirectional_iterator<base_base_t>
    {
        basic_iterator cpy{*this};
        --(*this);
        return cpy;
    }

    //!\brief Advances the iterator by skip positions.
    constexpr basic_iterator & operator+=(difference_type const skip)
        noexcept(noexcept(std::declval<base_t &>() += skip))
        requires std::random_access_iterator<base_base_t>
    {
        base_t::operator+=(skip);
        pos += skip;
        return *this;
    }

    //!\brief Advances the iterator by -skip positions.
    constexpr basic_iterator & operator-=(difference_type const skip)
        noexcept(noexcept(std::declval<base_t &>() -= skip))
        requires std::random_access_iterator<base_base_t>
    {
        base_t::operator-=(skip);
        pos -= skip;
        return *this;
    }
    //!\}

    /*!\name Comparison operators
     * \brief We define comparison against self and against the sentinel.
     * \{
     */

    //!\brief Checks whether `*this` is equal to `rhs`.
    constexpr bool operator==(basic_iterator const & rhs) const
        noexcept(!or_throw && noexcept(std::declval<base_base_t &>() == std::declval<base_base_t &>()))
        requires std::forward_iterator<base_base_t>
    {
        return this->base() == rhs.base();
    }

    //!\copydoc operator==()
    constexpr bool operator==(sentinel_type const & rhs) const
        noexcept(!or_throw && noexcept(std::declval<base_base_t const &>() == std::declval<sentinel_type const &>()))
    {
        if (pos >= max_pos)
            return true;

        if (this->base() == rhs)
        {
            if constexpr (or_throw)
                throw unexpected_end_of_input{"Reached end of input before designated size."};

            return true;
        }
        else
        {
            return false;
        }
    }

    //!\brief Checks whether `lhs` is equal to `rhs`.
    constexpr friend bool operator==(sentinel_type const & lhs, basic_iterator const & rhs)
        noexcept(noexcept(rhs == lhs))
    {
        return rhs == lhs;
    }

    //!\brief Checks whether `*this` is not equal to `rhs`.
    constexpr bool operator!=(sentinel_type const & rhs) const
        noexcept(noexcept(std::declval<basic_iterator &>() == rhs))
    {
        return !(*this == rhs);
    }

    //!\copydoc operator!=()
    constexpr bool operator!=(basic_iterator const & rhs) const
        noexcept(noexcept(std::declval<basic_iterator &>() == rhs))
        requires std::forward_iterator<base_base_t>
    {
        return !(*this == rhs);
    }

    //!\brief Checks whether `lhs` is not equal to `rhs`.
    constexpr friend bool operator!=(sentinel_type const & lhs, basic_iterator const & rhs)
        noexcept(noexcept(rhs != lhs))
    {
        return rhs != lhs;
    }
    //!\}

    /*!\name Reference/Dereference operators
     * \brief seqan3::detail::inherited_iterator_base operators are used unless specialised here.
     * \{
     */

    /*!\brief Accesses an element by index.
     * \param n Position relative to current location.
     * \return A reference to the element at relative location.
     */
    constexpr reference operator[](std::make_unsigned_t<difference_type> const n) const
        noexcept(noexcept(std::declval<base_base_t &>()[0]))
        requires std::random_access_iterator<base_base_t>
    {
        return base_t::operator[](n);
    }
    //!\}
};

// ============================================================================
//  take_fn (adaptor definition)
// ============================================================================

/*!\brief View adaptor definition for views::take and views::take_or_throw.
 * \tparam or_throw Whether to throw an exception when the input is exhausted before the end is reached.
 */
template <bool or_throw>
struct take_exactly_fn
{
    //!\brief Store the arguments and return a range adaptor closure object.
    constexpr auto operator()(size_t const size) const
    {
        return adaptor_from_functor{*this, size};
    }

    /*!\brief Type erase if possible and return view_take_exactly if not.
     * \returns An instance of std::span, std::basic_string_view, std::ranges::subrange or seqan3::detail::view_take_exactly.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, size_t target_size) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
                      "The views::take adaptor can only be passed viewable_ranges, i.e. Views or &-to-non-View.");

        // safeguard against wrong size
        if constexpr (std::ranges::sized_range<urng_t>)
        {
            if constexpr (or_throw)
            {
                if (target_size > std::ranges::size(urange))
                {
                    throw std::invalid_argument{"You are trying to construct a detail::take_exactly_or_throw from a "
                                                "range that is strictly smaller."};
                }
            }
            else
            {
                target_size = std::min<size_t>(target_size, std::ranges::size(urange));
            }
        }

        // string_view
        if constexpr (is_type_specialisation_of_v<std::remove_cvref_t<urng_t>, std::basic_string_view>)
        {
            // in standard
            return urange.substr(0, target_size);
        }
        // string const &
        else if constexpr (is_type_specialisation_of_v<std::remove_cvref_t<urng_t>, std::basic_string>
                           && std::is_const_v<std::remove_reference_t<urng_t>>)
        {
            // not in standard
            // seqan3::views::type_reduce does this too
            return std::basic_string_view{std::ranges::data(urange), target_size};
        }
        // contiguous
        else if constexpr (std::ranges::borrowed_range<urng_t> && std::ranges::contiguous_range<urng_t>
                           && std::ranges::sized_range<urng_t>)
        {
            // not in standard (special case for std::span in standard)
            // seqan3::views::type_reduce does this too
            return std::span{std::ranges::data(urange), target_size};
        }
        // random_access
        else if constexpr (std::ranges::borrowed_range<urng_t> && std::ranges::random_access_range<urng_t>
                           && std::ranges::sized_range<urng_t>)
        {
            // not in standard
            // seqan3::views::type_reduce does this too
            return std::ranges::subrange<std::ranges::iterator_t<urng_t>, std::ranges::iterator_t<urng_t>>{
                std::ranges::begin(urange),
                std::ranges::begin(urange) + target_size,
                target_size};
        }
        // our type
        else
        {
            return view_take_exactly<std::views::all_t<urng_t>, or_throw>{std::forward<urng_t>(urange), target_size};
        }
    }
};

} // namespace seqan3::detail

// ============================================================================
//  detail::take_exactly (adaptor instance definition)
// ============================================================================

namespace seqan3::detail
{
/*!\brief               A view adaptor that returns the first `size` elements from the underlying range (or less if the
 *                      underlying range is shorter); also provides size information.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \param[in] size      The target size of the view.
 * \returns             Up to `size` elements of the underlying range.
 * \ingroup io_views
 *
 * \details
 *
 * \header_file{seqan3/io/views/detail/take_exactly_view.hpp}
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |----------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                                        |
 * | std::ranges::forward_range       |                                       | *preserved*                                        |
 * | std::ranges::bidirectional_range |                                       | *preserved*                                        |
 * | std::ranges::random_access_range |                                       | *preserved*                                        |
 * | std::ranges::contiguous_range    |                                       | *preserved*                                        |
 * |                                  |                                       |                                                    |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                                       |
 * | std::ranges::view                |                                       | *guaranteed*                                       |
 * | std::ranges::sized_range         |                                       | ***guaranteed***                                   |
 * | std::ranges::common_range        |                                       | *preserved*                                        |
 * | std::ranges::output_range        |                                       | *preserved* except if `urng_t` is std::basic_string|
 * | std::ranges::borrowed_range      |                                       | *preserved*                                        |
 * | seqan3::const_iterable_range     |                                       | *preserved*                                        |
 * |                                  |                                       |                                                    |
 * | std::ranges::range_reference_t   |                                       | std::ranges::range_reference_t<urng_t>             |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * The difference to seqan3::views::take is that this view always exposes size information – even if the underlying
 * range is not sized. You should only use this if you know that the underlying range will always be
 * at least `size` long.
 *
 * For seqan3::detail::take_exactly if the underlying range is shorter than `size`, the behaviour is undefined.
 * seqan3::detail::take_exactly_or_throw is a safer alternative, because it throws an exception when an iterator before
 * the `size`-th one compares equal to the end sentinel; and it also throws on construction if it knows that the
 * underlying range is smaller.
 *
 * ### Example
 *
 * \include test/snippet/io/views/detail/take_exactly.cpp
 *
 * \hideinitializer
 */
inline constexpr auto take_exactly = take_exactly_fn<false>{};

// ============================================================================
//  detail::take_exactly_or_throw (adaptor instance definition)
// ============================================================================

/*!\brief A view adaptor that returns the first `size` elements from the underlying range and also exposes size
 *        information; throws if the underlying range is smaller than `size`.
 * \throws seqan3::unexpected_end_of_input If the underlying range is smaller than `size`.
 * \ingroup io_views
 *
 * \copydetails seqan3::detail::take_exactly
 * \hideinitializer
 */
inline constexpr auto take_exactly_or_throw = take_exactly_fn<true>{};
} // namespace seqan3::detail

namespace std::ranges
{
//!\cond
template <std::ranges::view urng_t, bool or_throw>
inline constexpr bool enable_borrowed_range<seqan3::detail::view_take_exactly<urng_t, or_throw>> =
    enable_borrowed_range<urng_t>;
//!\endcond
} // namespace std::ranges
