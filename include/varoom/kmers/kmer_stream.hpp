#ifndef VAROOM_KMERS_KMER_STREAM_HPP
#define VAROOM_KMERS_KMER_STREAM_HPP

#ifndef VAROOM_KMERS_NON_UNIQUE_KMER_STREAM_HPP
#include "varoom/kmers/non_unique_kmer_stream.hpp"
#endif

#ifndef VAROOM_KMERS_UNIQUE_KMER_STREAM_HPP
#include "varoom/kmers/unique_kmer_stream.hpp"
#endif

namespace varoom
{
    typedef std::vector<uint8_t> bytes;

    class basic_stream_writer
    {
    public:
        basic_stream_writer(bytes& p_unique, bytes& p_non_unique)
            : m_unique(p_unique), m_non_unique(p_non_unique)
        {
        }

        basic_stream_writer& push_back(const std::pair<uint64_t,uint64_t>& p_item)
        {
            if (p_item.second == 1)
            {
                m_unique.push_back(p_item.first);
            }
            else
            {
                m_non_unique.push_back(p_item);
            }
            return *this;
        }

    private:
        unique_kmer_stream_writer m_unique;
        non_unique_kmer_stream_writer m_non_unique;
    };

    class kmer_stream
    {
    public:
        template <typename LhsIterator, typename RhsIterator, typename Consumer>
        static void merge(LhsIterator& p_lhs, RhsIterator& p_rhs, Consumer& p_cons)
        {
            while (p_lhs.more() & p_rhs.more())
            {
                if ((*p_lhs).first < (*p_rhs).first)
                {
                    p_cons.push_back(*p_lhs);
                    ++p_lhs;
                    continue;
                }
                if ((*p_lhs).first > (*p_rhs).first)
                {
                    p_cons.push_back(*p_rhs);
                    ++p_rhs;
                    continue;
                }
                // kmers ==
                std::pair<uint64_t,uint64_t> itm = *p_lhs;
                itm.second += (*p_rhs).second;
                p_cons.push_back(itm);
                ++p_lhs;
                ++p_rhs;
            }
            while (p_lhs.more())
            {
                p_cons.push_back(*p_lhs);
                ++p_lhs;
            }
            while (p_rhs.more())
            {
                p_cons.push_back(*p_rhs);
                ++p_rhs;
            }
        }
    };
}
// namespace varoom

#endif // VAROOM_KMERS_KMER_STREAM_HPP
