#ifndef VAROOM_HGVS_HGVSG_HANDLER_HPP
#define VAROOM_HGVS_HGVSG_HANDLER_HPP

#include <cstdint>
#include <string>

namespace varoom
{
    namespace hgvs
    {
        class hgvsg_handler
        {
        public:
            virtual void sub(const std::string& p_chr, const std::int64_t p_pos, const std::string& p_ref, const std::string& p_alt) = 0;

            virtual void ins(const std::string& p_chr, const std::int64_t p_after, const std::int64_t p_before, const std::string& p_alt) = 0;

            virtual void del(const std::string& p_chr, const std::int64_t p_first, const std::int64_t p_last) = 0;

            virtual void delins(const std::string& p_chr, const std::int64_t p_first, const std::int64_t p_last, const std::string& p_seq) = 0;

            virtual void dup(const std::string& p_chr, const std::int64_t p_first, const std::int64_t p_last) = 0;

            virtual void inv(const std::string& p_chr, const std::int64_t p_first, const std::int64_t p_last) = 0;

        };
    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_HGVSG_HANDLER_HPP
