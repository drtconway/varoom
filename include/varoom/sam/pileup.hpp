#ifndef VAROOM_SAM_PILEUP_HPP
#define VAROOM_SAM_PILEUP_HPP

#ifndef VAROOM_SAM_HPP
#include "varoom/sam.hpp"
#endif

#include <map>

namespace varoom
{
    class sam_pileup
    {
    private:
        typedef std::map<std::uint32_t,std::map<std::string,size_t>> pileup;

    public:
        sam_pileup()
            : m_chr(""), m_pos(0)
        {
        }

        void add_alignment(const std::string& p_chr, const std::uint32_t& p_pos,
                           const std::string& p_seq, const std::string& p_cigar)
        {
            if (p_pos < m_pos)
            {
                throw std::domain_error("sam file not in sorted order.");
            }

            m_ops.clear();
            decode_cigar(p_cigar, m_ops);

            flush_upto(p_chr, p_pos);

            m_chr = p_chr;
            m_pos = p_pos;

            uint32_t i = 0;
            uint32_t j = 0;
            for (auto itr = m_ops.begin(); itr != m_ops.end(); ++itr)
            {
                switch (itr->first)
                {
                    case 'M':
                    {
                        for (uint32_t k = 0; k < itr->second; ++k)
                        {
                            char q = p_seq[i];
                            uint32_t p = p_pos + j;
                            m_pileup[p][atom(q)] += 1;
                            ++i;
                            ++j;
                        }
                        break;
                    }
                    case 'I':
                    {
                        uint32_t n = itr->second;
                        std::string q = std::string(p_seq.begin() + i, p_seq.begin() + i + n);
                        std::uint32_t p = p_pos + j;
                        m_pileup[p][q] += 1;
                        i += n;
                        break;
                    }
                    case 'D':
                    {
                        for (uint32_t k = 0; k < itr->second; ++k)
                        {
                            char q = '*';
                            std::uint32_t p = p_pos + j;
                            m_pileup[p][atom(q)] += 1;
                            ++j;
                        }
                        break;
                    }
                    case 'S':
                    {
                        i += itr->second;
                        j += itr->second;
                        break;
                    }
                    case 'H':
                    {
                        j += itr->second;
                        break;
                    }
                    default:
                    {
                        throw std::runtime_error("cigar op not implemented");
                    }
                }
            }
        }

        void end()
        {
            flush_upto("", 0);
        }

        virtual void output_pileup(const std::string& p_chr, const std::uint32_t& p_pos, const std::string& p_base, const size_t& p_count) = 0;
    private:
        void flush_upto(const std::string& p_chr, const std::uint32_t& p_pos)
        {
            while (true)
            {
                if (m_pileup.begin() == m_pileup.end() ||
                    (p_chr == m_chr && m_pileup.begin()->first >= p_pos))
                {
                    return;
                }
                const std::map<std::string,size_t>& counts = m_pileup.begin()->second;
                for (auto itr = counts.begin(); itr != counts.end(); ++itr)
                {
                    output_pileup(m_chr, m_pileup.begin()->first, itr->first, itr->second);
                }
                m_pileup.erase(m_pileup.begin());
            }
        }

        static const std::string& atom(char p_x)
        {
            static std::map<char,std::string> X;
            if (X.count(p_x) == 0)
            {
                X[p_x].push_back(p_x);
            }
            return X[p_x];
        }

        std::string m_chr;
        std::uint32_t m_pos;
        pileup m_pileup;
        std::vector<cigar_op> m_ops;
    };
}
// namespace varoom

#endif // VAROOM_SAM_PILEUP_HPP
