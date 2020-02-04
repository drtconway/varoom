#ifndef VAROOM_SEQ_GENOME_MAP_HPP
#define VAROOM_SEQ_GENOME_MAP_HPP

#include <string>
#include <unordered_map>
#include <boost/flyweight.hpp>
#include <nlohmann/json.hpp>

namespace varoom
{
    namespace seq
    {
        class genome_map
        {
        public:
            using chrom = boost::flyweight<std::string>;

            genome_map(const nlohmann::json& p_src)
            {
                size_t t = 0;
                for (auto itr = p_src.begin(); itr != p_src.end(); ++itr)
                {
                    chrom ch(itr.key());
                    size_t z = itr.value();
                    m_sizes[ch] = z;
                    m_toc[ch] = t;
                    m_index.push_back(ch);
                    m_starts.push_back(t);
                    t += z;
                }
                m_size = t;
            }

            size_t size() const
            {
                return m_size;
            }

            size_t size(const chrom& p_chr) const
            {
                return m_sizes.at(p_chr);
            }

            size_t chrom2genome(const chrom& p_chr, size_t p_pos) const
            {
                if (p_pos >= size(p_chr))
                {
                    throw std::runtime_error("position out of range for chromosome.");
                }
                return m_toc.at(p_chr) + p_pos;
            }

            std::pair<chrom,size_t> genome2chrom(size_t p_pos) const
            {
                if (p_pos >= size())
                {
                    throw std::runtime_error("position out of range");
                }
                auto itr = std::lower_bound(m_starts.begin(), m_starts.end(), p_pos);
                size_t j = itr - m_starts.begin();
                if (j > 0 && p_pos < *itr)
                {
                    j -= 1;
                }
                size_t pos = p_pos - m_starts[j];
                return std::make_pair(m_index[j], pos);
            }

        private:
            size_t m_size;
            std::unordered_map<chrom,size_t> m_sizes;
            std::unordered_map<chrom,size_t> m_toc;
            std::vector<chrom> m_index;
            std::vector<size_t> m_starts;
        };
    }
    // namespace seq
}
// namespace varoom

#endif // VAROOM_SEQ_GENOME_MAP_HPP
