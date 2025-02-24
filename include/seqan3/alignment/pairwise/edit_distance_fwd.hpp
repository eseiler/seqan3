// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Forwards for seqan3::edit_distance_unbanded related types.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
*/

#pragma once

#include <concepts>
#include <ranges>

#include <seqan3/alignment/configuration/align_config_min_score.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/detail/deferred_crtp_base.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/utility/detail/bits_of.hpp>

namespace seqan3::detail
{
template <typename word_t, typename score_t, bool is_semi_global, bool use_max_errors>
class edit_distance_score_matrix_full; //forward declaration

template <typename word_t, bool is_semi_global, bool use_max_errors>
class edit_distance_trace_matrix_full; //forward declaration

//!\brief Store no state for state_t.
template <typename state_t, typename...>
struct empty_state
{};

//!\brief If enabled is true state_t will be added to state_type.
template <bool enabled, typename state_t>
using enable_state_t = std::conditional_t<enabled, state_t, empty_state<state_t>>;

/*!\brief The same as std::conditional but for template template parameters.
 * \details
 * If `B` is true, <tt>selector<B, T, F>::template select</tt> inherits `T`, otherwise `F`.
 */
template <bool B, template <typename...> typename T, template <typename...> typename F>
struct selector
{
    //!\brief Depending on `B`, `select` is the template template parameter `T` or `F`.
    template <typename... args_t>
    struct select : public std::conditional_t<B, T<args_t...>, F<args_t...>>
    {};
};

/*!\brief The default traits type for the edit distance algorithm.
 * \ingroup alignment_pairwise
 */
template <std::ranges::viewable_range database_t,
          std::ranges::viewable_range query_t,
          typename align_config_t,
          typename is_semi_global_t,
          typename word_t = uint_fast64_t>
struct default_edit_distance_trait_type
{
    //!\brief The type of the alignment config.
    using align_config_type = std::remove_reference_t<align_config_t>;
    //!\brief The alignment algorithm traits over the alignment configuration type.
    using alignment_traits_type = alignment_configuration_traits<align_config_type>;
    //!\brief The type of one machine word.
    using word_type = word_t;
    static_assert(std::is_unsigned_v<word_type>, "the word type of edit_distance_unbanded must be unsigned.");
    static_assert(alignment_traits_type::has_output_configuration, "We assume the result type was configured.");
    //!\brief The type of the score.
    using score_type = typename alignment_traits_type::original_score_type;
    //!\brief The type of the database sequence.
    using database_type = std::remove_reference_t<database_t>;
    //!\brief The type of the query sequence.
    using query_type = std::remove_reference_t<query_t>;

    //!\brief The size of one machine word.
    static constexpr uint8_t word_size = bits_of<word_type>;
    static_assert(bits_of<word_type> <= 64u, "we assume at most uint64_t as word_type");

    //!\brief The type of an iterator of the database sequence.
    using database_iterator = std::ranges::iterator_t<database_type>;
    //!\brief The alphabet type of the query sequence.
    using query_alphabet_type = std::remove_reference_t<std::ranges::range_reference_t<query_type>>;
    //!\brief The alignment result type generated by the algorithm.
    using alignment_result_type = typename alignment_traits_type::alignment_result_type;
    //!\brief The alignment result value type.
    using result_value_type = typename alignment_result_value_type_accessor<alignment_result_type>::type;

    //!\brief When true the computation will use the ukkonen trick with the last active cell and bounds the error to
    //!       config.max_errors.
    static constexpr bool use_max_errors = align_config_type::template exists<align_cfg::min_score>();
    //!\brief Whether the alignment is a semi-global alignment or not.
    static constexpr bool is_semi_global = is_semi_global_t::value;
    //!\brief Whether the alignment is a global alignment or not.
    static constexpr bool is_global = !is_semi_global;
    //!\brief Whether the alignment configuration indicates to compute and/or store the score.
    static constexpr bool compute_score = true;
    //!\brief Whether the alignment configuration indicates to compute and/or store the alignment of the sequences.
    static constexpr bool compute_sequence_alignment = alignment_traits_type::compute_sequence_alignment;
    //!\brief Whether the alignment configuration indicates to compute and/or store the begin positions.
    static constexpr bool compute_begin_positions =
        alignment_traits_type::compute_begin_positions || compute_sequence_alignment;
    //!\brief Whether the alignment configuration indicates to compute and/or store the end positions.
    static constexpr bool compute_end_positions =
        alignment_traits_type::compute_end_positions || compute_begin_positions;
    //!\brief Whether the alignment configuration indicates to compute and/or store the score matrix.
    static constexpr bool compute_score_matrix = false;
    //!\brief Whether the alignment configuration indicates to compute and/or store the trace matrix.
    static constexpr bool compute_trace_matrix = compute_begin_positions || compute_sequence_alignment;
    //!\brief Whether the alignment configuration indicates to compute and/or store the score or trace matrix.
    static constexpr bool compute_matrix = compute_score_matrix || compute_trace_matrix;

    //!\brief The type of the trace matrix.
    using trace_matrix_type = edit_distance_trace_matrix_full<word_type, is_semi_global, use_max_errors>;
    //!\brief The type of the score matrix.
    using score_matrix_type = edit_distance_score_matrix_full<word_type, score_type, is_semi_global, use_max_errors>;
};

//!\brief A base class for edit_distance_unbanded.
template <bool enable_policy, template <typename...> typename policy_t, typename edit_traits, typename derived_t>
using edit_distance_base = invoke_deferred_crtp_base<
    deferred_crtp_base<selector<enable_policy, policy_t, empty_state>::template select, edit_traits>,
    derived_t>;

//!\cond
template <std::ranges::viewable_range database_t,
          std::ranges::viewable_range query_t,
          typename align_config_t,
          typename traits_t>
class edit_distance_unbanded; //forward declaration
//!\endcond

} // namespace seqan3::detail
