// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@.cpp.in

//![main]
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::rna4 letter1 = 'C'_dna4; // implicitly converted
    seqan3::rna4 letter2{};
    letter2 = 'C'_dna4; // implicitly converted
}
//![main]
