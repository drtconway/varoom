#ifndef VAROOM_HGVS_HGVSG_TO_HGVSC_HPP
#define VAROOM_HGVS_HGVSG_TO_HGVSC_HPP

#ifndef VAROOM_HGVS_HGVSG_HANDLER_HPP
#include "varoom/hgvs/hgvsg_handler.hpp"
#endif

#ifndef VAROOM_HGVS_HGVSC_HANDLER_HPP
#include "varoom/hgvs/hgvsc_handler.hpp"
#endif

#ifndef VAROOM_HGVS_HGVSC_LOCUS_HPP
#include "varoom/hgvs/hgvsc_locus.hpp"
#endif

#ifndef VAROOM_HGVS_TRANSCRIPT_HPP
#include "varoom/hgvs/transcript.hpp"
#endif

namespace varoom
{
    namespace hgvs
    {
        class hgvsg_to_hgvsc : public hgvsg_handler
        {
        public
            hgvsg_to_hgvsc(hgvsc_handler& p_handler)
                : m_handler(p_handler)
            {
            }
            
            virtual void sub(const std::string& p_chr, const std::int64_t p_pos, const std::string& p_ref, const std::string& p_alt)
            {
            }

            virtual void ins(const std::string& p_chr, const std::int64_t p_after, const std::int64_t p_before, const std::string& p_alt)
            {
            }

            virtual void del(const std::string& p_chr, const std::int64_t p_first, const std::int64_t p_last)
            {
            }

            virtual void delins(const std::string& p_chr, const std::int64_t p_first, const std::int64_t p_last, const std::string& p_seq)
            {
            }

            virtual void dup(const std::string& p_chr, const std::int64_t p_first, const std::int64_t p_last)
            {
            }

            virtual void inv(const std::string& p_chr, const std::int64_t p_first, const std::int64_t p_last)
            {
            }

        private:
            hgvsc_handler& m_handler;
        };
    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_HGVSG_TO_HGVSC_HPP
