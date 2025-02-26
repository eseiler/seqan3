// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Stream concepts.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <iosfwd>
#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3
{
/*!\interface seqan3::output_stream_over <>
 * \ingroup io_stream
 * \brief Concept for output streams.
 *
 * An object is an output stream if it inherits from the [std::ios_base](https://en.cppreference.com/w/cpp/io/ios_base)
 * and supports the (un)formatted output function (`operator<<`) for a l-value of a given `value_type`.
 * It further needs to define the public member types as described in the [STD](https://en.cppreference.com/w/cpp/io/basic_ios).
 */
//!\cond
template <typename stream_type, typename value_type>
concept output_stream_over =
    std::is_base_of_v<std::ios_base, std::remove_reference_t<stream_type>>
    && requires (stream_type & os, value_type & val) {
           typename std::remove_reference_t<stream_type>::char_type;
           typename std::remove_reference_t<stream_type>::traits_type;
           typename std::remove_reference_t<stream_type>::int_type;
           typename std::remove_reference_t<stream_type>::pos_type;
           typename std::remove_reference_t<stream_type>::off_type;

           {
               os << val
           } -> std::same_as<std::basic_ostream<typename std::remove_reference_t<stream_type>::char_type,
                                                typename std::remove_reference_t<stream_type>::traits_type> &>;
       };

template <typename stream_type>
concept output_stream = requires { typename std::remove_reference_t<stream_type>::char_type; }
                     && output_stream_over<stream_type, typename std::remove_reference_t<stream_type>::char_type>;
//!\endcond

/*!\name Requirements for seqan3::output_stream_over
 * \relates seqan3::output_stream_over
 * \brief You can expect these member types and the free function on all types that satisfy seqan3::output_stream_over.
 * \{
 */
/*!\fn      std::basic_ostream<char_type, traits_type> & operator<<(value_type val);
 * \brief   (un)-formatted output operator for the respective type on the underlying stream.
 * \param   val The value to write into the stream.
 * \returns A reference to a std::basic_ostream<char_type, traits_type>.
 *
 * \details
 *
 * The `char_type` and `traits_type` are inferred from the given `ostream`.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

/*!\typedef typename stream::char_type char_type
  * \brief Declares the associated char type.
  */

/*!\typedef typename stream::traits_type traits_type
  * \brief Declares the associated traits type.
  */

/*!\typedef typename stream::int_type int_type
  * \brief Declares the associated int type.
  */

/*!\typedef typename stream::pos_type pos_type
  * \brief Declares the associated pos type.
  */

/*!\typedef typename stream::off_type off_type
  * \brief Declares the associated off type.
  */
//!\}

/*!\interface seqan3::input_stream_over <>
 * \ingroup io_stream
 * \brief Concept for input streams.
 *
 * An object is an input stream if it inherits from the [std::ios_base](https://en.cppreference.com/w/cpp/io/ios_base)
 * and supports the (un)formatted input function (`operator>>`) for a l-value of a given `value_type`.
 * It further needs to define the public member types as described in the [STD](https://en.cppreference.com/w/cpp/io/basic_ios).
 */
//!\cond
template <typename stream_type, typename value_type>
concept input_stream_over =
    std::is_base_of_v<std::ios_base, std::remove_reference_t<stream_type>>
    && requires (stream_type & is, value_type & val) {
           typename std::remove_reference_t<stream_type>::char_type;
           typename std::remove_reference_t<stream_type>::traits_type;
           typename std::remove_reference_t<stream_type>::int_type;
           typename std::remove_reference_t<stream_type>::pos_type;
           typename std::remove_reference_t<stream_type>::off_type;

           {
               is >> val
           } -> std::same_as<std::basic_istream<typename std::remove_reference_t<stream_type>::char_type,
                                                typename std::remove_reference_t<stream_type>::traits_type> &>;
       };

template <typename stream_type>
concept input_stream = requires { typename std::remove_reference_t<stream_type>::char_type; }
                    && input_stream_over<stream_type, typename std::remove_reference_t<stream_type>::char_type>;
//!\endcond

/*!\name Requirements for seqan3::input_stream_over
 * \relates seqan3::input_stream_over
 * \brief You can expect these member types and the free function on all types that satisfy seqan3::input_stream_over.
 * \{
 */
/*!\fn      std::basic_istream<char_type, traits_type> & operator>>(value_type val);
 * \brief   (un)-formatted input operator for the respective type on the underlying stream.
 * \param   val The value to read from the stream.
 * \returns A reference to a std::basic_istream<char_type, traits_type>.
 *
 * \details
 *
 * The `char_type` and `traits_type` are inferred from the given `istream`.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

/*!\typedef typename stream::char_type char_type
  * \brief Declares the associated char type.
  */

/*!\typedef typename stream::traits_type traits_type
  * \brief Declares the associated traits type.
  */

/*!\typedef typename stream::int_type int_type
  * \brief Declares the associated int type.
  */

/*!\typedef typename stream::pos_type pos_type
  * \brief Declares the associated pos type.
  */

/*!\typedef typename stream::off_type off_type
  * \brief Declares the associated off type.
  */
//!\}

} // namespace seqan3
