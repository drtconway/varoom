#ifndef VAROOM_SEQ_FASTQ_FILTER_HPP
#define VAROOM_SEQ_FASTQ_FILTER_HPP

namespace varoom
{
    namespace seq
    {
        template <typename FastqStream>
        class fastq_filter
        {
        public:
            typedef FastqStream::item_type item_type;

            fastq_filter(FastqStream& p_src, std::function<bool(const fastq_read&)> p_pred)
                : m_src(p_src), m_pred(p_pred)
            {
            }

            bool more() const
            {
                return m_src.more();
            }

            const item_type& operator*() const
            {
                return *m_src;
            }

            fastq_filter& operator++()
            {
                next();
                return *this;
            }

        private:
            void next()
            {
                ++m_src;
                while (m_src.more() && !m_pred(*m_src))
                {
                    ++m_src;
                }
            }

            FastqStream& m_src;
            std::function<bool(const fastq_read&)> m_pred;
        };
    }
    // namespace seq
}
// namespace varoom

#endif // VAROOM_SEQ_FASTQ_FILTER_HPP
