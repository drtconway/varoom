#ifndef VAROOM_SEQ_LOCUS_STREAM_HPP
#define VAROOM_SEQ_LOCUS_STREAM_HPP

#include <string>
#include <boost/flyweight.hpp>

namespace varoom
{
    namespace seq
    {
        class locus_id
        {
        public:
            locus_id()
                : m_chr("<unspecified>"), m_ord(-1), m_pos(0)
            {
            }

            locus_id(const std::string& p_chr, uint32_t p_pos)
                : m_chr(p_chr), m_ord(-1), m_pos(p_pos)
            {
                auto itr = ord().find(p_chr);
                if (itr != ord().end())
                {
                    m_ord = itr->second;
                }
            }

            bool operator<(const locus_id& p_other) const
            {
                if (m_ord >= 0 && p_other.m_ord >= 0)
                {
                    return m_ord < p_other.m_ord || (m_ord == p_other.m_ord && m_pos < p_other.m_pos);
                }
                if (m_ord >= 0)
                {
                    return true;
                }
                if (p_other.m_ord >= 0)
                {
                    return false;
                }
                return m_chr < p_other.m_chr || (m_chr == p_other.m_chr && m_pos < p_other.m_pos);
            }

            bool operator==(const locus_id& p_other) const
            {
                if (m_ord >= 0 && p_other.m_ord >= 0)
                {
                    return m_ord == p_other.m_ord && m_pos == p_other.m_pos;
                }
                return m_chr == p_other.m_chr && m_pos == p_other.m_pos;
            }

            const std::string& chr() const
            {
                return m_chr;
            }

            const uint32_t pos() const
            {
                return m_pos;
            }

        private:
            static const std::map<std::string,int>& ord()
            {
                static std::map<std::string,int> m = {
                    {"1", 1}, {"chr1", 1}, {"2", 2}, {"chr2", 2},
                    {"3", 3}, {"chr3", 3}, {"4", 4}, {"chr4", 4},
                    {"5", 5}, {"chr5", 5}, {"6", 6}, {"chr6", 6},
                    {"7", 7}, {"chr7", 7}, {"8", 8}, {"chr8", 8},
                    {"9", 9}, {"chr9", 9}, {"10", 10}, {"chr10", 10},
                    {"11", 11}, {"chr11", 11}, {"12", 12}, {"chr12", 12},
                    {"13", 13}, {"chr13", 13}, {"14", 14}, {"chr14", 14},
                    {"15", 15}, {"chr15", 15}, {"16", 16}, {"chr16", 16},
                    {"17", 17}, {"chr17", 17}, {"18", 18}, {"chr18", 18},
                    {"19", 19}, {"chr19", 19}, {"20", 20}, {"chr20", 20},
                    {"21", 21}, {"chr21", 21}, {"22", 22}, {"chr22", 22},
                    {"X", 23}, {"chrX", 23}, {"Y", 24}, {"chrY", 24}
                };
                return m;
            }

            boost::flyweight<std::string> m_chr;
            int m_ord;
            uint32_t m_pos;
        };

        template <typename T>
        class locus_stream
        {
        public:
            locus_stream(const T& p_data)
                : m_more(true), m_data(p_data)
            {
            }

            bool more() const
            {
                return m_more;
            }

            const locus_id& locus() const
            {
                return m_locus;
            }

            const T& data() const
            {
                return m_data;
            }

            void operator++()
            {
                locus_id prev = m_locus;
                next();
                if (more() && !(prev < m_locus))
                {
                    throw std::runtime_error("mis-ordered loci");
                }
            }

        protected:
            virtual void next() = 0;

            bool m_more;
            locus_id m_locus;
            T m_data;
        };

        template <typename T,typename U,typename V>
        static void merge(T& p_lhs, U& p_rhs, V& p_dst)
        {
            while (p_lhs.more() && p_rhs.more())
            {
                const locus_id& l = p_lhs.locus();
                const locus_id& r = p_rhs.locus();
                if (l < r)
                {
                    p_dst.push_back(*p_lhs);
                    ++p_lhs;
                    continue;
                }
                if (r < l)
                {
                    p_dst.push_back(*p_rhs);
                    ++p_rhs;
                    continue;
                }
                p_dst.push_back(*p_lhs, *p_rhs);
                ++p_lhs;
                ++p_rhs;
            }
        }
    }
    // namespace seq
}
// namespace varoom

#endif // VAROOM_SEQ_LOCUS_STREAM_HPP
