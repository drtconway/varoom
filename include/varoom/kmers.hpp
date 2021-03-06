#ifndef VAROOM_KMERS_HPP
#define VAROOM_KMERS_HPP

#include <algorithm>
#include <bitset>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace varoom
{
    typedef std::uint64_t kmer;
    typedef std::pair<kmer,size_t> kmer_and_pos;

    namespace detail
    {
        inline size_t popcnt(const kmer& p_x)
        {
            return std::bitset<64>(p_x).count();
        }

        inline kmer rev_pairs(const kmer& p_x)
        {
            static constexpr kmer m2 = 0x3333333333333333ULL;
            static constexpr kmer m3 = 0x0F0F0F0F0F0F0F0FULL;
            static constexpr kmer m4 = 0x00FF00FF00FF00FFULL;
            static constexpr kmer m5 = 0x0000FFFF0000FFFFULL;
            static constexpr kmer m6 = 0x00000000FFFFFFFFULL;

            kmer x = p_x;
            x = ((x >> 2) & m2) | ((x & m2) << 2);
            x = ((x >> 4) & m3) | ((x & m3) << 4);
            x = ((x >> 8) & m4) | ((x & m4) << 8);
            x = ((x >> 16) & m5) | ((x & m5) << 16);
            x = ((x >> 32) & m6) | ((x & m6) << 32);
            return x;
        }
    }
    // namespace detail

    class kmers
    {
    public:
        static size_t ham(const kmer& p_x, const kmer& p_y)
        {
            static constexpr kmer m = 0x5555555555555555ULL;
            kmer z = p_x ^ p_y;
            return detail::popcnt((z | (z >> 1)) & m);
        }

        static kmer reverse_complement(const size_t& p_k, const kmer& p_x)
        {
            return detail::rev_pairs(~p_x) >> (64 - 2*p_k);
        }

        static kmer canonical(const size_t& p_k, const kmer& p_x)
        {
            kmer y = reverse_complement(p_k, p_x);
            return std::min(p_x, y);
        }

        static bool to_base(const char& p_c, kmer& p_x)
        {
            switch (p_c)
            {
                case 'A':
                case 'a':
                    p_x = 0;
                    return true;
                case 'C':
                case 'c':
                    p_x = 1;
                    return true;
                case 'G':
                case 'g':
                    p_x = 2;
                    return true;
                case 'T':
                case 't':
                    p_x = 3;
                    return true;
                default:
                    return false;
            }
        }

        static char to_char(const kmer& p_x)
        {
            return "ACGT"[p_x&3];
        }

        static void make(const std::string& p_seq, const size_t& p_k, std::function<void(kmer)> p_consumer)
        {
            const kmer M = (1ULL << (2*p_k)) - 1;
            kmer x = 0;
            size_t j = 0;
            for (auto s = p_seq.begin(); s != p_seq.end(); ++s)
            {
                kmer b;
                if (!to_base(*s, b))
                {
                    x = 0;
                    j = 0;
                    continue;
                }
                x = (x << 2) | b;
                ++j;
                if (j == p_k)
                {
                    p_consumer(x & M);
                    --j;
                }
            }
        }

        static void make(const std::string& p_seq, const size_t& p_k, std::vector<kmer>& p_kmers)
        {
            const kmer M = (1ULL << (2*p_k)) - 1;
            p_kmers.clear();
            kmer x = 0;
            size_t j = 0;
            for (auto s = p_seq.begin(); s != p_seq.end(); ++s)
            {
                kmer b;
                if (!to_base(*s, b))
                {
                    x = 0;
                    j = 0;
                    continue;
                }
                x = (x << 2) | b;
                ++j;
                if (j == p_k)
                {
                    p_kmers.push_back(x & M);
                    --j;
                }
            }
        }

        static void make(const std::string& p_seq, const size_t& p_k,
                         std::vector<kmer>& p_fwd,
                         std::vector<kmer>& p_rev)
        {
            const kmer M = (1ULL << (2*p_k)) - 1;
            const size_t S = 2*(p_k - 1);
            p_fwd.clear();
            p_rev.clear();
            kmer x = 0;
            kmer xb = 0;
            size_t j = 0;
            for (auto s = p_seq.begin(); s != p_seq.end(); ++s)
            {
                kmer b;
                if (!to_base(*s, b))
                {
                    x = 0;
                    j = 0;
                    continue;
                }
                kmer bb = 3-b;
                x = (x << 2) | b;
                xb = (xb >> 2) | (bb << S);
                ++j;
                if (j == p_k)
                {
                    p_fwd.push_back(x & M);
                    p_rev.push_back(xb);
                    --j;
                }
            }
        }

        template <class X>
        static void make(const std::string& p_seq, const size_t& p_k, X p_acceptor)
        {
            static_assert(std::is_convertible<X, std::function<void(kmer_and_pos)>>::value,
                          "applicative requires a function type void(kmer_and_pos)");
            const kmer M = (1ULL << (2*p_k)) - 1;
            kmer x = 0;
            size_t i = 0;
            size_t j = 0;
            for (auto s = p_seq.begin(); s != p_seq.end(); ++s, ++i)
            {
                kmer b;
                if (!to_base(*s, b))
                {
                    x = 0;
                    j = 0;
                    continue;
                }
                x = (x << 2) | b;
                ++j;
                if (j == p_k)
                {
                    p_acceptor(kmer_and_pos(x & M, i - p_k + 1));
                    --j;
                }
            }
        }

        static void make(const std::string& p_seq, const size_t& p_k,
                         std::vector<kmer_and_pos>& p_fwd,
                         std::vector<kmer_and_pos>& p_rev)
        {
            const kmer M = (1ULL << (2*p_k)) - 1;
            const size_t S = 2*(p_k - 1);
            p_fwd.clear();
            p_rev.clear();
            kmer x = 0;
            kmer xb = 0;
            size_t i = 0;
            size_t ib = p_seq.size() - 1;
            size_t j = 0;
            for (auto s = p_seq.begin(); s != p_seq.end(); ++s, ++i, --ib)
            {
                kmer b;
                if (!to_base(*s, b))
                {
                    x = 0;
                    j = 0;
                    continue;
                }
                kmer bb = 3-b;
                x = (x << 2) | b;
                xb = (xb >> 2) | (bb << S);
                ++j;
                if (j == p_k)
                {
                    p_fwd.push_back(kmer_and_pos(x & M, i - p_k + 1));
                    p_rev.push_back(kmer_and_pos(xb, ib));
                    --j;
                }
            }
        }

        static void render(const size_t& p_k, const kmer& p_x, std::string& p_res)
        {
            p_res.clear();
            kmer x = p_x;
            for (size_t i = 0; i < p_k; ++i)
            {
                p_res.push_back(to_char(x));
                x >>= 2;
            }
            std::reverse(p_res.begin(), p_res.end());
        }

        static std::string render(const size_t& p_k, const kmer& p_x)
        {
            std::string s;
            render(p_k, p_x, s);
            return s;
        }
    };
    // namespace kmers
}
// namespace varoom

#endif // VAROOM_KMERS_HPP
