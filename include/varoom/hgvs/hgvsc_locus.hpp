#ifndef VAROOM_HGVS_HGVSC_LOCUS_HPP
#define VAROOM_HGVS_HGVSC_LOCUS_HPP

#ifndef VAROOM_HGVS_HGVS_POSITIONS_HPP
#include "varoom/hgvs/hgvs_positions.hpp"
#endif

namespace varoom
{
    namespace hgvs
    {
        enum tx_zone { UTR5, CODING, INTRON, UTR3 };

        struct hgvsc_locus
        {
            tx_zone zone;
            hgvsc_position pos;
            hgvsc_relative_position rel;

            hgvsc_locus()
                : zone(CODING), pos(1), rel(0)
            {
            }

            hgvsc_locus(tx_zone p_zone, const hgvsc_position& p_pos, const hgvsc_relative_position& p_rel)
                : zone(p_zone), pos(p_pos), rel(p_rel)
            {
            }
        };
    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_HGVSC_LOCUS_HPP
