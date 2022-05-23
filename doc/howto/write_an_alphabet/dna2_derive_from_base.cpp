// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

//! [dna2]
#include <array> // std::array

#include <seqan3/alphabet/alphabet_base.hpp>            // alphabet_base
#include <seqan3/alphabet/concept.hpp>                  // alphabet concept checks
#include <seqan3/utility/char_operations/transform.hpp> // seqan3::to_lower

// derive from alphabet_base
struct dna2 : public seqan3::alphabet_base<dna2, 2>
{
private:
    // make the base class a friend so it can access the tables:
    friend seqan3::alphabet_base<dna2, 2>;

    // map 0 => 'S' and 1 => 'W'
    static constexpr char_type rank_to_char(rank_type const rank)
    {
        // via a lookup table
        return rank_to_char_table[rank];
        // or via an arithmetic expression
        return rank == 1 ? 'W' : 'S';
    }

    static constexpr rank_type char_to_rank(char_type const chr)
    {
        // via a lookup table
        using index_t = std::make_unsigned_t<char_type>;
        return char_to_rank_table[static_cast<index_t>(chr)];
        // or via an arithmetic expression
        return seqan3::to_lower(chr) == 'w' ? 1 : 0;
    }

private:
    // === lookup-table implementation detail ===

    // map 0 => 'S' and 1 => 'W'
    static constexpr char_type rank_to_char_table[alphabet_size]{'S', 'W'};

    static constexpr std::array<rank_type, 256> char_to_rank_table{
        // initialise with an immediately evaluated lambda expression:
        []() constexpr {std::array<rank_type, 256> ret{}; // initialise all values with 0 (=> 'S')
    ret['W'] = 1;                                         // only 'W' and 'w' result in rank 1
    ret['w'] = 1;
    return ret;
}()
}
;
}
;

// check the concepts
static_assert(seqan3::alphabet<dna2>);          // ok
static_assert(seqan3::writable_alphabet<dna2>); // ok
//! [dna2]
