#ifndef VAROOM_HGVS_HGVSG_FORMATTER_HPP
#define VAROOM_HGVS_HGVSG_FORMATTER_HPP

#ifndef VAROOM_HGVS_HGVSG_HANDLER_HPP
#include "varoom/hgvs/hgvsg_handler.hpp"
#endif

#include <functional>

namespace varoom
{
    namespace hgvs
    {
        class hgvsg_formatter : public hgvsg_handler
        {
        public:
            hgvsg_formatter(std::function<void(const std::string&)> p_consumer)
                : m_consumer(p_consumer)
            {
            }

            virtual void sub(const std::string& p_chr, const std::int64_t p_pos, const std::string& p_ref, const std::string& p_alt)
            {
                m_consumer(p_chr + ":g." + std::to_string(p_pos) + p_ref + ">" + p_alt);
            }

            virtual void ins(const std::string& p_chr, const std::int64_t p_after, const std::int64_t p_before, const std::string& p_alt)
            {
                m_consumer(p_chr + ":g." + std::to_string(p_after) + "_" + std::to_string(p_before) + "ins" + p_alt);
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
            std::function<void(const std::string&)> m_consumer;
        };
    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_HGVSG_FORMATTER_HPP

