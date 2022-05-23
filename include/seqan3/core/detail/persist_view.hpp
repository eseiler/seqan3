// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::detail::persist.
 */

#pragma once

#include <algorithm>
#include <concepts>
#include <seqan3/std/ranges>

#include <seqan3/core/detail/iterator_traits.hpp>
#include <seqan3/core/range/detail/adaptor_base.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/utility/range/concept.hpp>
#include <seqan3/utility/type_traits/detail/transformation_trait_or.hpp>

namespace seqan3::detail
{

// ============================================================================
//  view_persist
// ============================================================================

/*!\brief The type returned by seqan3::detail::persist.
 * \tparam urng_t The type of the underlying range, must model std::ranges::input_range.
 * \implements std::ranges::view
 * \implements std::ranges::random_access_range
 * \implements std::ranges::sized_range
 * \ingroup core
 *
 * \details
 *
 * Note that most members of this class are generated by ranges::view_interface which is not yet documented here.
 */
template <std::ranges::input_range urng_t>
class view_persist : public std::ranges::view_interface<view_persist<urng_t>>
{
private:
    //!\brief Shared storage of the underlying range.
    std::shared_ptr<urng_t> urange;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    view_persist() noexcept = default;                                               //!< Defaulted.
    constexpr view_persist(view_persist const & rhs) noexcept = default;             //!< Defaulted.
    constexpr view_persist(view_persist && rhs) noexcept = default;                  //!< Defaulted.
    constexpr view_persist & operator=(view_persist const & rhs) noexcept = default; //!< Defaulted.
    constexpr view_persist & operator=(view_persist && rhs) noexcept = default;      //!< Defaulted.
    ~view_persist() noexcept = default;                                              //!< Defaulted.

    /*!\brief Construct from another range.
     * \param[in] _urange The underlying range.
     */
    view_persist(urng_t && _urange) : urange{new urng_t{std::move(_urange)}}
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
    auto begin() noexcept
    {
        return std::ranges::begin(*urange);
    }

    //!\copydoc begin()
    auto begin() const noexcept
        requires const_iterable_range<urng_t>
    {
        return std::ranges::cbegin(*urange);
    }

    /*!\brief Returns an iterator to the element following the last element of the range.
     * \returns Iterator or sentinel to the end.
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
    auto end() noexcept
    {
        return std::ranges::end(*urange);
    }

    //!\copydoc end()
    auto end() const noexcept
        requires const_iterable_range<urng_t>
    {
        return std::ranges::cend(*urange);
    }
    //!\}
};

//!\brief Template argument type deduction guide that strips references.
//!\relates seqan3::detail::view_persist
template <typename urng_t>
view_persist(urng_t &&) -> view_persist<std::remove_reference_t<urng_t>>;

// ============================================================================
//  persist_fn (adaptor definition)
// ============================================================================

//![adaptor_def]
/*!\brief View adaptor definition for detail::persist.
 *!\ingroup core
 */
class persist_fn : public adaptor_base<persist_fn>
{
private:
    //!\brief Type of the CRTP-base.
    using base_t = adaptor_base<persist_fn>;

public:
    //!\brief Inherit the base class's Constructors.
    using base_t::base_t;

private:
    //!\brief Befriend the base class so it can call impl().
    friend base_t;

    /*!\brief       For ranges that are viewable, delegate to std::views::all.
     * \returns     An instance of std::views::all.
     */
    template <std::ranges::viewable_range urng_t>
    static auto impl(urng_t && urange)
    {
        return std::views::all(std::forward<urng_t>(urange));
    }

    /*!\brief       For ranges that are not views and not lvalue-references, call view_persist's constructor.
     * \returns     An instance of seqan3::detail::view_persist.
     */
    template <std::ranges::range urng_t>
    static auto impl(urng_t && urange)
    {
        static_assert(!std::is_lvalue_reference_v<urng_t>, "BUG: lvalue-reference in persist_fn::impl().");
        return view_persist{std::move(urange)};
    }
};
//![adaptor_def]

// ============================================================================
//  detail::persist (adaptor instance definition)
// ============================================================================

/*!\brief               A view adaptor that wraps rvalue references of non-views.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range wrapped in a view (even if it doesn't model std::ranges::viewable_range).
 * \ingroup core
 *
 * \details
 *
 * \header_file{seqan3/core/detail/persist_view.hpp}
 *
 * For ranges that model std::ranges::viewable_range, this adaptor just returns std::views::all. However this adaptor
 * can also take ranges that are not "viewable", e.g. temporaries of containers. It wraps them in a shared pointer
 * internally so all view requirements like constant copy are satisfied. However construction and copying might be
 * slightly slower, because of reference counting.
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
 * | std::ranges::viewable_range      | ***not required***                    | *guaranteed*                                       |
 * | std::ranges::view                |                                       | *guaranteed*                                       |
 * | std::ranges::sized_range         |                                       | *preserved*                                        |
 * | std::ranges::common_range        |                                       | *preserved*                                        |
 * | std::ranges::output_range        |                                       | *preserved*                                        |
 * | seqan3::const_iterable_range     |                                       | *preserved*                                        |
 * |                                  |                                       |                                                    |
 * | std::ranges::range_reference_t   |                                       | std::ranges::range_reference_t<urng_t>             |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * \include test/snippet/core/detail/persist_view.cpp
 *
 * \hideinitializer
 */
inline constexpr auto persist = persist_fn{};

} // namespace seqan3::detail
