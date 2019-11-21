#ifndef VAROOM_HGVS_HGVSC_FORMATTER_HPP
#define VAROOM_HGVS_HGVSC_FORMATTER_HPP

#ifndef VAROOM_HGVS_HGVSC_HANDLER_HPP
#include "varoom/hgvs/hgvsc_handler.hpp"
#endif

#include <functional>
#include <boost/format.hpp>

namespace varoom
{
    namespace hgvs
    {
        class hgvsc_formatter : public hgvsc_handler
        {
        public:
            hgvsc_formatter(std::function<void(const std::string&)> p_consumer)
                : m_consumer(p_consumer)
            {
            }

            virtual void sub(const std::string& p_chr, const hgvsc_locus& p_pos, const std::string& p_ref, const std::string& p_alt)
            {
                m_consumer(boost::str(boost::format("%s:c.%s%s>%s") % p_chr % p_pos.str() % p_ref % p_alt));
            }

            virtual void ins(const std::string& p_chr, const hgvsc_locus& p_after, const hgvsc_locus& p_before, const std::string& p_alt)
            {
                m_consumer(boost::str(boost::format("%s:c.%s_%sins%s") % p_chr % p_after.str() % p_before.str() % p_alt));
            }

            virtual void del(const std::string& p_chr, const hgvsc_locus& p_first, const hgvsc_locus& p_last)
            {
                if (p_first == p_last)
                {
                    m_consumer(boost::str(boost::format("%s:c.%sdel") % p_chr % p_first.str()));
                }
                else
                {
                    m_consumer(boost::str(boost::format("%s:c.%s_%sdel") % p_chr % p_first.str() % p_last.str()));
                }
            }

            virtual void delins(const std::string& p_chr, const hgvsc_locus& p_first, const hgvsc_locus& p_last, const std::string& p_seq)
            {
                if (p_first == p_last)
                {
                    m_consumer(boost::str(boost::format("%s:c.%sdelins%s") % p_chr % p_first.str() % p_seq));
                }
                else
                {
                    m_consumer(boost::str(boost::format("%s:c.%s_%sdelins%s") % p_chr % p_first.str() % p_last.str() % p_seq));
                }
            }

            virtual void dup(const std::string& p_chr, const hgvsc_locus& p_first, const hgvsc_locus& p_last)
            {
                if (p_first == p_last)
                {
                    m_consumer(boost::str(boost::format("%s:c.%sdup") % p_chr % p_first.str()));
                }
                else
                {
                    m_consumer(boost::str(boost::format("%s:c.%s_%sdup") % p_chr % p_first.str() % p_last.str()));
                }
            }

            virtual void inv(const std::string& p_chr, const hgvsc_locus& p_first, const hgvsc_locus& p_last)
            {
                m_consumer(boost::str(boost::format("%s:c.%s_%sinv") % p_chr % p_first.str() % p_last.str()));
            }

        private:
            std::function<void(const std::string&)> m_consumer;
        };
    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_HGVSC_FORMATTER_HPP

