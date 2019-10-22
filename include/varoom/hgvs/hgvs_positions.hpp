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
        struct gpos_0_tag;
        typedef scalar<gpos_0_tag> gpos_0;

        struct gpos_1_tag;
        typedef scalar<gpos_1_tag> gpos_1;

        struct txpos_0_tag;
        typedef scalar<txpos_0_tag> txpos_0;

        struct txpos_1_tag;
        typedef scalar<txpos_1_tag> txpos_1;
    }
    // namespace hgvs

    template<>
    struct scalar_conversions<hgvs::gpos_0_tag,hgvs::gpos_1_tag>
    {
        static scalar<hgvs::gpos_1_tag> cast(const scalar<hgvs::gpos_0_tag>& p_x)
        {
            return scalar<hgvs::gpos_1_tag>(p_x() + 1);
        }
    };

    template<>
    struct scalar_conversions<hgvs::gpos_1_tag,hgvs::gpos_0_tag>
    {
        static scalar<hgvs::gpos_0_tag> cast(const scalar<hgvs::gpos_1_tag>& p_x)
        {
            return scalar<hgvs::gpos_0_tag>(p_x() - 1);
        }
    };

    template<>
    struct scalar_conversions<hgvs::txpos_0_tag,hgvs::txpos_1_tag>
    {
        static scalar<hgvs::txpos_1_tag> cast(const scalar<hgvs::txpos_0_tag>& p_x)
        {
            if (p_x() >= 0)
            {
                return scalar<hgvs::txpos_1_tag>(p_x() + 1);
            }
            else
            {
                return scalar<hgvs::txpos_1_tag>(p_x());
            }
        }
    };

    template<>
    struct scalar_conversions<hgvs::txpos_1_tag,hgvs::txpos_0_tag>
    {
        static scalar<hgvs::txpos_0_tag> cast(const scalar<hgvs::txpos_1_tag>& p_x)
        {
            if (p_x() > 0)
            {
                return scalar<hgvs::txpos_0_tag>(p_x() - 1);
            }
            if (p_x() < 0)
            {
                return scalar<hgvs::txpos_0_tag>(p_x());
            }
            throw std::invalid_argument("1-offset transcript positions cannot be 0");
        }
    };
}
// namespace varoom

#endif // VAROOM_HGVS_HGVS_POSITIONS_HPP
