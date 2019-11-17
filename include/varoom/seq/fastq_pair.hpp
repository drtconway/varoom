#ifndef VAROOM_SEQ_FASTQ_PAIR_HPP
#define VAROOM_SEQ_FASTQ_PAIR_HPP

#ifndef VAROOM_SEQ_FASTQ_HPP
#include "varoom/seq/fastq.hpp"
#endif

#include <functional>

namespace varoom
{
    namespace seq
    {
        template <typename LhsStream, typename RhsStream = LhsStream>
        class fastq_pair
        {
        public:
            typedef std::pair<fastq_reader::item_result_type,fastq_reader::item_result_type> item_type;
            typedef item_type item_result_type;

            fastq_pair(LhsStream& p_lhs, RhsStream& p_rhs)
                : m_lhs(p_lhs), m_rhs(p_rhs)
            {
            }

            bool more() const
            {
                return m_lhs.more();
            }

            item_result_type operator*() const
            {
                return item_type(*m_lhs, *m_rhs);
            }

            fastq_pair& operator++()
            {
                next();
                return *this;
            }

        private:
            void next()
            {
                ++m_lhs;
                ++m_rhs;
                if (m_lhs.more() != m_rhs.more())
                {
                    throw std::runtime_error("fastq streams not of equal length");
                }
            }

            LhsStream& m_lhs;
            RhsStream& m_rhs;
        };
    }
    // namespace seq
}
// namespace varoom

#endif // VAROOM_SEQ_FASTQ_PAIR_HPP
