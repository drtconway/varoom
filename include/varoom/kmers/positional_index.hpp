#ifndef VAROOM_KMERS_KMER_ACCUMULATOR_HPP
#define VAROOM_KMERS_KMER_ACCUMULATOR_HPP

#ifndef VAROOM_KMERS_HPP
#include "varoom/kmers.hpp"
#endif

#include <iostream>
#include <unordered_map>

namespace varoom
{
    class positional_index
    {
    public:
        typedef size_t seq_num;
        typedef int32_t seq_pos;
        typedef std::pair<seq_num,seq_pos> seq_num_and_pos;
        typedef std::pair<seq_pos,size_t> seq_pos_and_count;

        positional_index(const size_t& p_k)
            : m_k(p_k), m_num_seqs(0)
        {
        }

        seq_num add(const std::string& p_label, const std::string& p_seq, bool p_both_strands = false)
        {
            seq_num n = add(p_seq, p_both_strands);
            m_label_index[n] = p_label;
            return n;
        }

        seq_num add(const std::string& p_seq, bool p_both_strands = false)
        {
            std::vector<kmer_and_pos> fwd;
            std::vector<kmer_and_pos> rev;

            seq_num n = m_num_seqs++;
            kmers::make(p_seq, m_k, fwd, rev);

            add_kmers(n, fwd);
            if (p_both_strands)
            {
                add_kmers(n, rev);
            }
            return n;
        }

        bool hits(const std::vector<kmer_and_pos>& p_xs, std::unordered_map<seq_num, std::unordered_map<seq_pos,size_t>>& p_res) const
        {
            std::vector<seq_num_and_pos> hx;
            for (size_t i = 0; i < p_xs.size(); ++i)
            {
                kmer x = p_xs[i].first;
                int32_t p = p_xs[i].second;

                size_t j = hx.size();
                hits(x, hx);

                for (; j < hx.size(); ++j)
                {
                    hx[j].second -= p;
                }
            }
            if (hx.size() == 0)
            {
                return false;
            }

            std::sort(hx.begin(), hx.end());

            seq_num pn = 0;
            seq_pos pp = 0;
            size_t cx = 0;
            for (size_t i = 0; i < hx.size(); ++i)
            {
                seq_num n = hx[i].first;
                seq_pos p = hx[i].second;
                if (n != pn || p != pp)
                {
                    if (cx > 0)
                    {
                        p_res[pn][pp] += cx;
                    }
                    pn = n;
                    pp = p;
                    cx = 0;
                }
                cx += 1;
            }
            if (cx > 0)
            {
                p_res[pn][pp] += cx;
            }

            return true;
        }

        bool hits(const kmer& p_x, std::vector<seq_num_and_pos>& p_res) const
        {
            auto itr = m_index.find(p_x);
            if (itr == m_index.end())
            {
                return false;
            }
            p_res.insert(p_res.end(), itr->second.begin(), itr->second.end());
            return true;
        }

    private:
        static std::string hex(const kmer& p_x)
        {
            std::string r;
            uint64_t x = p_x;
            for (size_t i = 0; i < 64; i += 4)
            {
                r.push_back("0123456789ABCDEF"[x & 0xf]);
                x >>= 4;
            }
            return std::string(r.rbegin(), r.rend());
        }

        void validate_kmer(const kmer& p_x) const
        {
            uint64_t v = p_x >> (2*m_k);
            if (v != 0)
            {
                std::cout << "bad! " << kmers::render(m_k, p_x) << '\t' << hex(v) << std::endl;
                throw std::runtime_error("invalid kmer");
            }
        }

        void add_kmers(const size_t p_n, const std::vector<kmer_and_pos>& p_xs)
        {
            for (size_t i = 0; i < p_xs.size(); ++i)
            {
                kmer x = p_xs[i].first;
                size_t p = p_xs[i].second;
                m_index[x].push_back(seq_num_and_pos(p_n, p));
            }
        }

        const size_t m_k;
        size_t m_num_seqs;
        std::unordered_map<seq_num,std::string> m_label_index;
        std::unordered_map<kmer,std::vector<seq_num_and_pos>> m_index;
    };
}
// namespace varoom

#endif // VAROOM_KMERS_KMER_ACCUMULATOR_HPP
