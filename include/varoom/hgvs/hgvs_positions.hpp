#ifndef VAROOM_HGVS_HGVS_POSITIONS_HPP
#define VAROOM_HGVS_HGVS_POSITIONS_HPP

#ifndef VAROOM_UTIL_STRONG_TYPEDEF_HPP
#include "varoom/util/strong_typedef.hpp"
#endif

#include <stdexcept>

namespace varoom
{
    namespace hgvs
    {
        // The type genomic_locus is for referring to a 0-based
        // position in a genomic reference sequence.
        //
        struct genomic_locus : varoom::strong_typedef<genomic_locus, uint64_t>, varoom::integer_arithmetic<genomic_locus>
        {
            using strong_typedef::strong_typedef;
        };

        // The type hgvsg_locus is for referring to a 1-based
        // position in a genomic reference sequence, a la HGVSg.
        //
        struct hgvsg_locus : varoom::strong_typedef<hgvsg_locus, uint64_t>, varoom::integer_arithmetic<hgvsg_locus>
        {
            using strong_typedef::strong_typedef;
        };

        // The type tx_position is for referring to a 0-based position
        // in a transcript relative to the start codon for coding transcripts,
        // or the transcription start for non-coding transcripts.
        //
        struct tx_position : varoom::strong_typedef<tx_position, int64_t>, varoom::integer_arithmetic<tx_position>
        {
            using strong_typedef::strong_typedef;
        };

        // The type hgvsc_position is for referring to the 1-based position 
        // in a transcript relative to the start codon for coding transcripts,
        // or the transcription start for non-coding transcripts.
        // Note that negative positions refer to the 5' UTR, and there is no
        // directed representation for intronic positions. This is why there
        // is a separate type hgvsc_locus.
        //
        struct hgvsc_position : varoom::strong_typedef<hgvsc_position, int64_t>, varoom::integer_arithmetic<hgvsc_position>
        {
            using strong_typedef::strong_typedef;
        };
    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_HGVS_POSITIONS_HPP
