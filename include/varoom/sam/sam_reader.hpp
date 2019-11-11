#ifndef VAROOM_SAM_SAM_READER_HPP
#define VAROOM_SAM_SAM_READER_HPP

#include <istream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/iterator_range.hpp>

#ifndef VAROOM_SAM_SAM_ALIGNMENT_HPP
#include "varoom/sam/sam_alignment.hpp"
#endif

namespace varoom
{
    class sam_reader
    {
    public:
        sam_reader(std::istream& p_in)
            : m_in(p_in), m_more(true)
        {
            next();
        }

        bool more() const
        {
            return m_more;
        }

        const sam_alignment& operator*() const
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
            while (m_more)
            {
                if (!std::getline(m_in, m_line))
                {
                    m_more = false;
                    return;
                }
                if (starts_with(m_line, '@'))
                {
                    // Metadata
                    continue;
                }
                break;
            }

            boost::algorithm::trim(m_line);
            subtext st(m_line);
            st.split('\t', m_parts);
            m_curr.name.clear();
            m_curr.name.insert(m_curr.name.end(), m_parts[0].first, m_parts[0].second);
            m_curr.flags = boost::lexical_cast<uint32_t>(boost::make_iterator_range(m_parts[1].first, m_parts[1].second));
            m_curr.chr.clear();
            m_curr.chr.insert(m_curr.chr.end(), m_parts[2].first, m_parts[2].second);
            m_curr.pos = boost::lexical_cast<uint32_t>(boost::make_iterator_range(m_parts[3].first, m_parts[3].second));
            m_curr.mapq = boost::lexical_cast<uint32_t>(boost::make_iterator_range(m_parts[4].first, m_parts[4].second));
            m_curr.cigar.clear();
            m_curr.cigar.insert(m_curr.cigar.end(), m_parts[5].first, m_parts[5].second);
            m_curr.mate_chr.clear();
            m_curr.mate_chr.insert(m_curr.mate_chr.end(), m_parts[6].first, m_parts[6].second);
            m_curr.mate_pos = boost::lexical_cast<uint32_t>(boost::make_iterator_range(m_parts[7].first, m_parts[7].second));
            m_curr.tlen = boost::lexical_cast<uint32_t>(boost::make_iterator_range(m_parts[8].first, m_parts[8].second));
            m_curr.seq.clear();
            m_curr.seq.insert(m_curr.seq.end(), m_parts[9].first, m_parts[9].second);
            m_curr.qual.clear();
            m_curr.qual.insert(m_curr.qual.end(), m_parts[10].first, m_parts[10].second);
        }

        static bool starts_with(const std::string& p_str, char p_ch)
        {
            return p_str.size() > 0 && p_str.front() == p_ch;
        }

        std::istream& m_in;
        bool m_more;
        std::string m_line;
        std::vector<subtext> m_parts;
        sam_alignment m_curr;
    };

}
// namespace varoom

#endif // VAROOM_SAM_SAM_READER_HPP
