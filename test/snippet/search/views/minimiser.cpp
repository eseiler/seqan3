#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/syncmer.hpp>

using namespace seqan3::literals;

int main()
{
    static constexpr size_t const syncmer_t{0u};
    static constexpr size_t const syncmer_k{3u};
    static constexpr size_t const syncmer_s{2u};

    std::vector<seqan3::dna4> text{"GGCAAGTGACA"_dna4};
    auto text_begin = text.begin();

    // auto hashes = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{syncmer_s}});
    // seqan3::debug_stream << hashes << '\n';

    auto syncmer = text | seqan3::views::syncmer({.kmer_size = syncmer_k, .smer_size = syncmer_s, .offset = syncmer_t});

    for (auto syncmer_it = syncmer.begin(); syncmer_it != syncmer.end(); ++syncmer_it)
    {
        auto text_pos = syncmer_it.base() - text_begin - syncmer_k + 1;
        seqan3::debug_stream << "Text Pos: " << text_pos << " Text: " << std::span{text_begin + text_pos, syncmer_k}
                             << " Value: " << *syncmer_it << '\n';
    }
}
