#ifndef VAROOM_HGVS_HGVSC_LOCUS_HPP
#define VAROOM_HGVS_HGVSC_LOCUS_HPP

#ifndef VAROOM_HGVS_HGVS_POSITIONS_HPP
#include "varoom/hgvs/hgvs_positions.hpp"
#endif

#include <boost/format.hpp>

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

            bool operator==(const hgvsc_locus& p_other) const
            {
                return zone == p_other.zone && pos == p_other.pos && rel == p_other.rel;
            }

            std::string str() const
            {
                switch (zone)
                {
                    case UTR5:
                    {
                        if (static_cast<int64_t>(rel) < 0)
                        {
                            return boost::str(boost::format("%d%d") % static_cast<int64_t>(pos) % static_cast<int64_t>(rel));
                        }
                        else if (static_cast<int64_t>(rel) == 0)
                        {
                            return boost::str(boost::format("%d") % static_cast<int64_t>(pos));
                        }
                        else
                        {
                            return boost::str(boost::format("%d+%d") % static_cast<int64_t>(pos) % static_cast<int64_t>(rel));
                        }
                    }
                    case CODING:
                    {
                        return boost::str(boost::format("%d") % static_cast<int64_t>(pos));
                    }
                    case INTRON:
                    {
                        if (static_cast<int64_t>(rel) < 0)
                        {
                            return boost::str(boost::format("%d%d") % static_cast<int64_t>(pos) % static_cast<int64_t>(rel));
                        }
                        else
                        {
                            return boost::str(boost::format("%d+%d") % static_cast<int64_t>(pos) % static_cast<int64_t>(rel));
                        }
                    }
                    case UTR3:
                    {
                        if (static_cast<int64_t>(rel) < 0)
                        {
                            return boost::str(boost::format("*%d%d") % static_cast<int64_t>(pos) % static_cast<int64_t>(rel));
                        }
                        else if (static_cast<int64_t>(rel) == 0)
                        {
                            return boost::str(boost::format("*%d") % static_cast<int64_t>(pos));
                        }
                        else
                        {
                            return boost::str(boost::format("*%d+%d") % static_cast<int64_t>(pos) % static_cast<int64_t>(rel));
                        }
                    }
                }
                throw std::runtime_error("internal logic error");
            }
        };
    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_HGVSC_LOCUS_HPP
