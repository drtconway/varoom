#ifndef VAROOM_KMERS_NON_UNIQUE_KMER_STREAM_HPP
#define VAROOM_KMERS_NON_UNIQUE_KMER_STREAM_HPP

#ifndef VAROOM_UTIL_CODEC8_HPP
#include "varoom/util/codec8.hpp"
#endif

namespace varoom
{
    namespace kmers
    {
        typedef std::vector<uint8_t> bytes;

        class non_unique_kmer_stream_reader
        {
        public:
            non_unique_kmer_stream_reader(const bytes& p_bytes)
                : m_curr(p_bytes.begin()), m_end(p_bytes.end()), m_xac(0, 0), m_more(true)
            {
                next();
            }

            bool more() const
            {
                return m_more;
            }

            const std::pair<uint64_t,uint64_t>& operator*() const
            {
                return m_xac;
            }

            non_unique_kmer_stream_reader& operator++()
            {
                next();
                return *this;
            }

        private:
            void next()
            {
                if (m_curr == m_end)
                {
                    m_more = false;
                    return;
                }
                m_xac.first += varoom::util::codec8::decode(m_curr);
                m_xac.second = varoom::util::codec8::decode(m_curr);
            }

            bytes::const_iterator m_curr;
            bytes::const_iterator m_end;
            std::pair<uint64_t,uint64_t> m_xac;
            bool m_more;
        };

        class non_unique_kmer_stream_writer
        {
        public:
            non_unique_kmer_stream_writer(bytes& p_bytes)
                : m_bytes(p_bytes), m_px(0)
            {
            }

            non_unique_kmer_stream_writer& push_back(const std::pair<uint64_t,uint64_t>& p_itm)
            {
                return push_back(p_itm.first, p_itm.second);
            }

            non_unique_kmer_stream_writer& push_back(kmer p_x, size_t p_count)
            {
                varoom::util::codec8::encode(p_x - m_px, m_bytes);
                varoom::util::codec8::encode(p_count, m_bytes);
                m_px = p_x;
                return *this;
            }

        private:
            bytes& m_bytes;
            uint64_t m_px;
        };
    }
    // namespace kmers
}
// namespace varoom

#endif // VAROOM_KMERS_NON_UNIQUE_KMER_STREAM_HPP
