#ifndef VAROOM_HGVS_HGVS_POSITIONS_HPP
#define VAROOM_HGVS_HGVS_POSITIONS_HPP

#ifndef VAROOM_UTIL_SCALAR_HPP
#include "varoom/util/scalar.hpp"
#endif

#include <stdexcept>

namespace varoom
{
    namespace hgvs
    {
        // The type genomic_locus is for referring to a 0-based
        // position in a genomic reference sequence.
        //
        struct genomic_locus_tag;
        typedef scalar<genomic_locus_tag> genomic_locus;

        // The type hgvsg_locus is for referring to a 1-based
        // position in a genomic reference sequence, a la HGVSg.
        //
        struct hgvsg_locus_tag;
        typedef scalar<hgvsg_locus_tag> hgvsg_locus;

        // The type tx_position is for referring to a 0-based position
        // in a transcript relative to the start codon for coding transcripts,
        // or the transcription start for non-coding transcripts.
        //
        struct tx_position_tag;
        typedef scalar<tx_position_tag> tx_position;

        // The type hgvsc_position is for referring to the 1-based position 
        // in a transcript relative to the start codon for coding transcripts,
        // or the transcription start for non-coding transcripts.
        // Note that negative positions refer to the 5' UTR, and there is no
        // directed representation for intronic positions. This is why there
        // is a separate type hgvsc_locus.
        //
        struct hgvsc_position_tag;
        typedef scalar<hgvsc_position_tag> hgvsc_position;
    }
    // namespace hgvs

    template<>
    struct scalar_conversions<hgvs::genomic_locus_tag,hgvs::hgvsg_locus_tag>
    {
        static scalar<hgvs::hgvsg_locus_tag> cast(const scalar<hgvs::genomic_locus_tag>& p_x)
        {
            return scalar<hgvs::hgvsg_locus_tag>(p_x() + 1);
        }
    };

    template<>
    struct scalar_conversions<hgvs::hgvsg_locus_tag,hgvs::genomic_locus_tag>
    {
        static scalar<hgvs::genomic_locus_tag> cast(const scalar<hgvs::hgvsg_locus_tag>& p_x)
        {
            return scalar<hgvs::genomic_locus_tag>(p_x() - 1);
        }
    };

    template<>
    struct scalar_conversions<hgvs::tx_position_tag,hgvs::hgvsc_position_tag>
    {
        static scalar<hgvs::hgvsc_position_tag> cast(const scalar<hgvs::tx_position_tag>& p_x)
        {
            if (p_x() >= 0)
            {
                return scalar<hgvs::hgvsc_position_tag>(p_x() + 1);
            }
            else
            {
                return scalar<hgvs::hgvsc_position_tag>(p_x());
            }
        }
    };

    template<>
    struct scalar_conversions<hgvs::hgvsc_position_tag,hgvs::tx_position_tag>
    {
        static scalar<hgvs::tx_position_tag> cast(const scalar<hgvs::hgvsc_position_tag>& p_x)
        {
            if (p_x() > 0)
            {
                return scalar<hgvs::tx_position_tag>(p_x() - 1);
            }
            if (p_x() < 0)
            {
                return scalar<hgvs::tx_position_tag>(p_x());
            }
            throw std::invalid_argument("1-offset transcript positions cannot be 0");
        }
    };
}
// namespace varoom

#endif // VAROOM_HGVS_HGVS_POSITIONS_HPP
