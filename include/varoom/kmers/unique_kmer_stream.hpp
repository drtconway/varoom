#ifndef VAROOM_KMERS_UNIQUE_KMER_STREAM_HPP
#define VAROOM_KMERS_UNIQUE_KMER_STREAM_HPP

#ifndef VAROOM_UTIL_CODEC8_HPP
#include "varoom/util/codec8.hpp"
#endif

namespace varoom
{
    typedef std::vector<uint8_t> bytes;

    class unique_kmer_stream_reader
    {
    public:
        unique_kmer_stream_reader(const bytes& p_bytes)
            : m_curr(p_bytes.begin()), m_end(p_bytes.end()), m_xac(0, 1), m_more(true)
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

        unique_kmer_stream_reader& operator++()
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
        }

        bytes::const_iterator m_curr;
        bytes::const_iterator m_end;
        std::pair<uint64_t,uint64_t> m_xac;
        bool m_more;
    };

    class unique_kmer_stream_writer
    {
    public:
        unique_kmer_stream_writer(bytes& p_bytes)
            : m_bytes(p_bytes), m_px(0)
        {
        }

        unique_kmer_stream_writer& push_back(kmer p_x)
        {
            varoom::util::codec8::encode(p_x - m_px, m_bytes);
            m_px = p_x;
            return *this;
        }

    private:
        bytes& m_bytes;
        uint64_t m_px;
    };
}
// namespace varoom

#endif // VAROOM_KMERS_UNIQUE_KMER_STREAM_HPP
