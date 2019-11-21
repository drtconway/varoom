#ifndef VAROOM_HGVS_HGVSC_HANDLER_HPP
#define VAROOM_HGVS_HGVSC_HANDLER_HPP

#include "varoom/hgvs/hgvsc_locus.hpp"

namespace varoom
{
    namespace hgvs
    {
        class hgvsc_handler
        {
        public:
            virtual void sub(const std::string& p_chr, const hgvsc_locus& p_pos, const std::string& p_ref, const std::string& p_alt) = 0;

            virtual void ins(const std::string& p_chr, const hgvsc_locus& p_after, const hgvsc_locus& p_before, const std::string& p_alt) = 0;

            virtual void del(const std::string& p_chr, const hgvsc_locus& p_first, const hgvsc_locus& p_last) = 0;

            virtual void delins(const std::string& p_chr, const hgvsc_locus& p_first, const hgvsc_locus& p_last, const std::string& p_seq) = 0;

            virtual void dup(const std::string& p_chr, const hgvsc_locus& p_first, const hgvsc_locus& p_last) = 0;

            virtual void inv(const std::string& p_chr, const hgvsc_locus& p_first, const hgvsc_locus& p_last) = 0;

        };
    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_HGVSC_HANDLER_HPP

