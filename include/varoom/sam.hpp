#ifndef VAROOM_SAM_HPP
#define VAROOM_SAM_HPP

#include <cstdint>
#include <istream>
#include <regex>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/iterator_range.hpp>

#ifndef VAROOM_UTIL_SUBTEXT_HPP
#include "varoom/util/subtext.hpp"
#endif

namespace varoom
{
    struct sam_flags
    {
        static constexpr std::uint32_t paired = 1;
        static constexpr std::uint32_t proper_pair = 2;
        static constexpr std::uint32_t unmapped = 4;
        static constexpr std::uint32_t mate_unmapped = 8;
        static constexpr std::uint32_t reverse = 16;
        static constexpr std::uint32_t mate_reverse = 32;
        static constexpr std::uint32_t read_1 = 64;
        static constexpr std::uint32_t read_2 = 128;
        static constexpr std::uint32_t secondary = 256;
        static constexpr std::uint32_t qcfail = 512;
        static constexpr std::uint32_t dup = 1024;
        static constexpr std::uint32_t supplementary = 2048;

        static bool is_paired(const std::uint32_t& p_flgs)
        {
            return p_flgs & paired;
        }

        static bool is_proper_pair(const std::uint32_t& p_flgs)
        {
            return p_flgs & proper_pair;
        }

        static bool is_unmapped(const std::uint32_t& p_flgs)
        {
            return p_flgs & unmapped;
        }

        static bool is_mate_unmapped(const std::uint32_t& p_flgs)
        {
            return p_flgs & mate_unmapped;
        }

        static bool is_reverse(const std::uint32_t& p_flgs)
        {
            return p_flgs & reverse;
        }

        static bool is_mate_reverse(const std::uint32_t& p_flgs)
        {
            return p_flgs & mate_reverse;
        }

        static bool is_read_1(const std::uint32_t& p_flgs)
        {
            return p_flgs & read_1;
        }

        static bool is_read_2(const std::uint32_t& p_flgs)
        {
            return p_flgs & read_2;
        }

        static bool is_secondary(const std::uint32_t& p_flgs)
        {
            return p_flgs & secondary;
        }

        static bool is_qcfail(const std::uint32_t& p_flgs)
        {
            return p_flgs & qcfail;
        }

        static bool is_dup(const std::uint32_t& p_flgs)
        {
            return p_flgs & dup;
        }

        static bool is_supplementary(const std::uint32_t& p_flgs)
        {
            return p_flgs & supplementary;
        }
    };

    struct sam_alignment
    {
        std::string name;
        std::uint32_t flags;
        std::string chr;
        std::uint32_t pos;
        std::uint32_t mapq;
        std::string cigar;
        std::string mate_chr;
        std::uint32_t mate_pos;
        std::uint32_t tlen;
        std::string seq;
        std::string qual;
    };

    typedef std::pair<char,std::uint32_t> cigar_op;

    void decode_cigar(const std::string& p_cigar, std::vector<cigar_op>& p_ops)
    {
        const std::regex cig_regex("([0-9]+)([DHIMNPSX=])");

        p_ops.clear();
        auto cig_begin = std::sregex_iterator(p_cigar.begin(), p_cigar.end(), cig_regex);
        auto cig_end = std::sregex_iterator();
        for (auto itr = cig_begin; itr != cig_end; ++itr)
        {
            std::smatch m = *itr;
            uint32_t n = boost::lexical_cast<uint32_t>(m[1]);
            char c = *(m[2].first);
            p_ops.push_back(cigar_op(c,n));
        }
    }

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
            m_parts.clear();
            subtext st(m_line);
            st.split('\t', m_parts);
            m_curr.name.clear();
            m_curr.name.insert(m_curr.name.end(), m_parts[0].first, m_parts[0].second);
            m_curr.flags = boost::lexical_cast<std::uint32_t>(boost::make_iterator_range(m_parts[1].first, m_parts[1].second));
            m_curr.chr.clear();
            m_curr.chr.insert(m_curr.chr.end(), m_parts[2].first, m_parts[2].second);
            m_curr.pos = boost::lexical_cast<std::uint32_t>(boost::make_iterator_range(m_parts[3].first, m_parts[3].second));
            m_curr.mapq = boost::lexical_cast<std::uint32_t>(boost::make_iterator_range(m_parts[4].first, m_parts[4].second));
            m_curr.cigar.clear();
            m_curr.cigar.insert(m_curr.cigar.end(), m_parts[5].first, m_parts[5].second);
            m_curr.mate_chr.clear();
            m_curr.mate_chr.insert(m_curr.mate_chr.end(), m_parts[6].first, m_parts[6].second);
            m_curr.mate_pos = boost::lexical_cast<std::uint32_t>(boost::make_iterator_range(m_parts[7].first, m_parts[7].second));
            m_curr.tlen = boost::lexical_cast<std::uint32_t>(boost::make_iterator_range(m_parts[8].first, m_parts[8].second));
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

#endif // VAROOM_SAM_HPP
