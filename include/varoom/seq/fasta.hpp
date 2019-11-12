#ifndef VAROOM_SEQ_FASTA_HPP
#define VAROOM_SEQ_FASTA_HPP

#include <istream>
#include <string>
#include <boost/algorithm/string/trim.hpp>

#ifndef VAROOM_UTIL_TEXT_HPP
#include "varoom/util/text.hpp"
#endif

namespace varoom
{
    namespace seq
    {
        typedef std::pair<std::string,std::string> fasta_read;

        class fasta_reader
        {
        public:
            fasta_reader(std::istream& p_in)
                : m_in(p_in), m_more(true)
            {
                next();
            }

            bool more() const
            {
                return m_more;
            }

            const fasta_read& operator*() const
            {
                return m_curr;
            }

            void operator++()
            {
                next();
            }

        private:
            void next()
            {
                if (m_next.size() == 0 && !std::getline(m_in, m_next))
                {
                    m_more = false;
                    return;
                }
                if (!varoom::text::starts_with(m_next, '>'))
                {
                    m_more = false;
                    throw std::domain_error("missing header line.");
                }

                m_curr.first.clear();
                m_curr.second.clear();

                boost::algorithm::trim(m_next);
                m_curr.first.insert(m_curr.first.end(), m_next.begin() + 1, m_next.end()); // drop the >

                while (std::getline(m_in, m_next) && !varoom::text::starts_with(m_next, '>'))
                {
                    // Technically, FASTA data may contain whitespace, so we should use code like:
                    //      str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
                    // but in practice this doesn't occur in cases we're interested in.
                    boost::algorithm::trim(m_next);
                    m_curr.second += m_next;
                }
            }

            std::istream& m_in;
            bool m_more;
            std::string m_next;
            fasta_read m_curr;
        };
    }
    // namespace seq
}
// namespace varoom

#endif // VAROOM_SEQ_FASTA_HPP
