// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_trace_matrix_full.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/matrix/detail/alignment_matrix_column_major_range_base.hpp>
#include <seqan3/alignment/matrix/detail/alignment_trace_matrix_base.hpp>
#include <seqan3/alignment/matrix/detail/alignment_trace_matrix_proxy.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief An alignment traceback matrix storing the entire traceback matrix.
 * \tparam trace_t          The type of the trace directions.
 * \tparam coordinate_only  A boolean flag indicating if only a seqan3::alignment_coordinate should be generated.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * This implementation allocates the full traceback matrix using quadratic memory. The matrix allows access to the
 * underlying values through a range based interface. An iterator over the traceback matrix iterates in
 * column-major-order over the traceback matrix. Dereferencing an iterator returns a view over the current matrix
 * column. The value type is a pair over a seqan3::alignment_coordinate and
 * the seqan3::detail::alignment_trace_matrix_proxy, which gives a unified access to the respective matrix cells
 * as needed by the standard alignment algorithm. The matrix is modelled as std::ranges::input_range since the
 * alignment algorithm iterates only once over the complete matrix to calculate the values.
 *
 * ### Only computing the coordinates
 *
 * Sometimes it is desired to only get access to the alignment coordinates. This can be achieved by setting
 * `coordinate_only = true`. In this case no memory will be allocated and only an internal state is maintained to
 * generate the alignment coordinates.
 */
template <typename trace_t, bool coordinate_only = false>
class alignment_trace_matrix_full :
    protected alignment_trace_matrix_base<trace_t>,
    public alignment_matrix_column_major_range_base<alignment_trace_matrix_full<trace_t, coordinate_only>>
{
private:
    static_assert(std::same_as<trace_t, trace_directions> || simd_concept<trace_t>,
                  "Value type must either be a trace_directions object or a simd vector.");

    //!\brief The base class for data storage.
    using matrix_base_t = alignment_trace_matrix_base<trace_t>;
    //!\brief The base class for iterating over the matrix.
    using range_base_t = alignment_matrix_column_major_range_base<alignment_trace_matrix_full<trace_t, coordinate_only>>;

    //!\brief Befriend the range base class.
    friend range_base_t;

protected:
    using typename matrix_base_t::element_type;
    using typename matrix_base_t::coordinate_type;
    using typename range_base_t::alignment_column_type;
    //!\brief The proxy type when accessing a traceback matrix cell.
    using proxy_type = std::conditional_t<coordinate_only,
                                          detail::ignore_t,
                                          alignment_trace_matrix_proxy<trace_t>>;

    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::column_data_view_type
    using column_data_view_type = std::conditional_t<coordinate_only,
                                        decltype(std::views::iota(coordinate_type{}, coordinate_type{})),
                                        decltype(std::views::zip(std::declval<std::span<element_type>>(),
                                                                std::declval<std::span<element_type>>(),
                                                                std::views::iota(coordinate_type{}, coordinate_type{})))>;

public:
    /*!\name Associated types
     * \{
     */
    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::value_type
    using value_type = std::pair<coordinate_type, proxy_type>;
    //!\brief Same as value type.
    using reference = value_type;
    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::iterator
    using iterator = typename range_base_t::iterator;
    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::sentinel
    using sentinel = typename range_base_t::sentinel;
    using typename matrix_base_t::size_type;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alignment_trace_matrix_full() = default; //!< Defaulted.
    constexpr alignment_trace_matrix_full(alignment_trace_matrix_full const &) = default; //!< Defaulted.
    constexpr alignment_trace_matrix_full(alignment_trace_matrix_full &&) = default; //!< Defaulted.
    constexpr alignment_trace_matrix_full & operator=(alignment_trace_matrix_full const &) = default; //!< Defaulted.
    constexpr alignment_trace_matrix_full & operator=(alignment_trace_matrix_full &&) = default; //!< Defaulted.
    ~alignment_trace_matrix_full() = default; //!< Defaulted.

    /*!\brief Construction from two ranges.
     * \tparam first_sequence_t  The first range type; must model std::ranges::forward_range.
     * \tparam second_sequence_t The second range type; must model std::ranges::forward_range.
     *
     * \param[in] first  The first range.
     * \param[in] second The second range.
     * \param[in] initial_value The value to initialise the matrix with. Default initialised if not specified.
     *
     * \details
     *
     * Obtains the sizes of the passed ranges in order to allocate the traceback matrix.
     * If `coordinate_only` is set to `true`, nothing will be allocated.
     */
    template <std::ranges::forward_range first_sequence_t, std::ranges::forward_range second_sequence_t>
    constexpr alignment_trace_matrix_full(first_sequence_t && first,
                                          second_sequence_t && second,
                                          [[maybe_unused]] trace_t const initial_value = trace_t{})
    {
        matrix_base_t::num_cols = static_cast<size_type>(std::ranges::distance(first) + 1);
        matrix_base_t::num_rows = static_cast<size_type>(std::ranges::distance(second) + 1);

        if constexpr (!coordinate_only)
        {
            matrix_base_t::data.resize(matrix_base_t::num_rows * matrix_base_t::num_cols);
            matrix_base_t::cache_left.resize(matrix_base_t::num_rows, initial_value);
        }
    }
    //!\}

private:
    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::initialise_column
    constexpr alignment_column_type initialise_column(size_type const column_index) noexcept
    {
        coordinate_type row_begin{column_index_type{column_index}, row_index_type{0u}};
        coordinate_type row_end{column_index_type{column_index}, row_index_type{matrix_base_t::num_rows}};
        if constexpr (coordinate_only)
        {
            return alignment_column_type{*this, column_data_view_type{std::views::iota(std::move(row_begin),
                                                                                      std::move(row_end))}};
        }
        else
        {
            size_type index = matrix_base_t::num_rows * column_index;
            auto col = std::views::zip(std::span<element_type>{std::addressof(matrix_base_t::data[index]),
                                                              matrix_base_t::num_rows},
                                      std::span<element_type>{matrix_base_t::cache_left},
                                      std::views::iota(std::move(row_begin), std::move(row_end)));
            return alignment_column_type{*this, column_data_view_type{col}};
        }
    }

    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::make_proxy
    template <std::random_access_iterator iter_t>
    constexpr value_type make_proxy(iter_t host_iter) noexcept
    {
        if constexpr (coordinate_only)
        {
            return {*host_iter, std::ignore};
        }
        else
        {
            return {std::get<2>(*host_iter),
                    proxy_type{std::get<0>(*host_iter),
                               std::get<1>(*host_iter),
                               std::get<1>(*host_iter),
                               matrix_base_t::cache_up}};
        }
    }
};

} // namespace seqan3::detail