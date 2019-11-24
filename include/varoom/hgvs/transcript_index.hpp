#ifndef VAROOM_HGVS_TRANSCRIPT_INDEX_HPP
#define VAROOM_HGVS_TRANSCRIPT_INDEX_HPP

#ifndef VAROOM_HGVS_TRANSCRIPT_HPP
#include "varoom/hgvs/transcript.hpp"
#endif

#include <unordered_map>

namespace varoom
{
    namespace hgvs
    {
        typedef std::shared_ptr<transcript> transcript_ptr;

        class transcript_index
        {
        public:
            transcript_index()
            {
            }

            void add(const transcript& p_tx)
            {
                m_index[p_tx.chr()].push_back(p_tx);
            }

            void overlapping(const std::string& p_chr, const genomic_locus& p_loc,
                             std::vector<std::reference_wrapper<const transcript>>& p_txs,
                             size_t p_window = 250) const
            {
                genomic_locus w(p_window);

                auto itr = m_index.find(p_chr);
                if (itr == m_index.end())
                {
                    return;
                }

                const std::vector<transcript>&  txs = itr->second;

                for (size_t i = 0; i < txs.size(); ++i)
                {
                    const transcript& tx = txs[i];
                    if (tx.tx_begin() - w <= p_loc && p_loc < tx.tx_end() + w)
                    {
                        p_txs.push_back(std::cref(tx));
                    }
                }
            }

        private:
            std::unordered_map<std::string,std::vector<transcript>> m_index;
        };
    }
    // namespace hgvs
}
// namespace varoom


#endif // VAROOM_HGVS_TRANSCRIPT_INDEX_HPP
