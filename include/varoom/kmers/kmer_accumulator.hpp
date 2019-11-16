#ifndef VAROOM_KMERS_KMER_ACCUMULATOR_HPP
#define VAROOM_KMERS_KMER_ACCUMULATOR_HPP

#ifndef VAROOM_KMERS_HPP
#include "varoom/kmers.hpp"
#endif

#ifndef VAROOM_KMERS_KMER_STREAM_HPP
#include "varoom/kmers/kmer_stream.hpp"
#endif

#include <deque>
#include <vector>

namespace varoom
{
    namespace kmers
    {
        class kmer_accumulator
        {
        public:
            kmer_accumulator()
            {
            }

            kmer_accumulator& push_back(const kmer& p_kmer)
            {
                m_buffer.push_back(p_kmer);
                return *this;
            }

            kmer_accumulator& operator+=(const std::vector<kmer>& p_kmers)
            {
                m_buffer.insert(m_buffer.end(), p_kmers.begin(), p_kmers.end());
                return *this;
            }

            template <typename Vis>
            void visit(Vis& p_vis)
            {
                flush();
                merge();
                if (m_unique.size() == 0 && m_non_unique.size() == 0)
                {
                    return;
                }
                if (m_unique.size() == 0)
                {
                    for (non_unique_kmer_stream_reader itr(m_non_unique.front()); itr.more(); ++itr)
                    {
                        p_vis.push_back(*itr);
                    }
                    return;
                }
                if (m_non_unique.size() == 0)
                {
                    for (unique_kmer_stream_reader itr(m_unique.front()); itr.more(); ++itr)
                    {
                        p_vis.push_back(*itr);
                    }
                    return;
                }
                unique_kmer_stream_reader xsItr(m_unique.front());
                non_unique_kmer_stream_reader ysItr(m_non_unique.front());
                varoom::kmers::merge(xsItr, ysItr, p_vis);
            }

        private:
            void flush()
            {
                bytes xs;
                bytes ys;
                {
                    basic_stream_writer dst(xs, ys);

                    uint64_t prev = 0;
                    uint64_t count = 0;
                    std::sort(m_buffer.begin(), m_buffer.end());
                    for (auto itr = m_buffer.begin(); itr != m_buffer.end(); ++itr)
                    {
                        uint64_t x = *itr;
                        if (x != prev)
                        {
                            if (count > 0)
                            {
                                std::pair<uint64_t,uint64_t> v(prev, count);
                                dst.push_back(v);
                            }
                            prev = x;
                            count = 1;
                        }
                        else
                        {
                            count += 1;
                        }
                    }
                    if (count > 0)
                    {
                        std::pair<uint64_t,uint64_t> v(prev, count);
                        dst.push_back(v);
                    }
                }
                if (xs.size() > 0)
                {
                    m_unique.push_back(xs);
                }
                if (ys.size() > 0)
                {
                    m_non_unique.push_back(ys);
                }
                m_buffer.clear();
            }

            void merge()
            {
                // Merge unique streams.
                //
                // Some k-mers will become non-unique
                //
                while (m_unique.size() > 1)
                {
                    bytes xs;
                    std::swap(xs, m_unique.front());
                    m_unique.pop_front();

                    bytes ys;
                    std::swap(ys, m_unique.front());
                    m_unique.pop_front();

                    unique_kmer_stream_reader xsItr(xs);
                    unique_kmer_stream_reader ysItr(ys);

                    merge(xsItr, ysItr);
                }

                // Merge non-unique streams.
                //
                // No k-mer can become unique!
                //
                while (m_non_unique.size() > 1)
                {
                    bytes xs;
                    std::swap(xs, m_non_unique.front());
                    m_non_unique.pop_front();

                    bytes ys;
                    std::swap(ys, m_non_unique.front());
                    m_non_unique.pop_front();

                    non_unique_kmer_stream_reader xsItr(xs);
                    non_unique_kmer_stream_reader ysItr(ys);

                    merge(xsItr, ysItr);
                }

                // Now merge the unique k-mers against the non-unique ones.
                //
                if (m_unique.size() > 0 && m_non_unique.size() > 0)
                {
                    bytes xs;
                    std::swap(xs, m_unique.front());
                    m_unique.pop_front();

                    bytes ys;
                    std::swap(ys, m_non_unique.front());
                    m_non_unique.pop_front();

                    unique_kmer_stream_reader xsItr(xs);
                    non_unique_kmer_stream_reader ysItr(ys);

                    merge(xsItr, ysItr);
                }
            }

            template <typename LhsIterator, typename RhsIterator>
            void merge(LhsIterator& p_lhs, RhsIterator& p_rhs)
            {
                bytes xs;
                bytes ys;
                {
                    basic_stream_writer dst(xs, ys);
                    varoom::kmers::merge(p_lhs, p_rhs, dst);
                }
                if (xs.size() > 0)
                {
                    m_unique.push_back(xs);
                }
                if (ys.size() > 0)
                {
                    m_non_unique.push_back(ys);
                }
            }

            std::vector<kmer> m_buffer;
            std::deque<bytes> m_unique;
            std::deque<bytes> m_non_unique;
        };
    }
    // namespace kmers
}
// namespace varoom


#endif // VAROOM_KMERS_KMER_ACCUMULATOR_HPP
