#ifndef VAROOM_XLR_HASH_HPP
#define VAROOM_XLR_HASH_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace varoom
{
    template<unsigned long long X, unsigned long long L, unsigned long long R>
    void hash_block(std::uint64_t& p_h)
    {
        p_h *= X;
        p_h ^= (p_h << L) | (p_h >> (64 - L));
        p_h ^= (p_h >> R) | (p_h << (64 - R));
    }

    template<unsigned long long X, unsigned long long L, unsigned long long R, unsigned long long Y, unsigned long long... More>
    void hash_block(std::uint64_t& p_h)
    {
        hash_block<X,L,R>(p_h);
        hash_block<Y,More...>(p_h);
    }

    void hash_block_d(std::uint64_t& p_h, const std::uint64_t& p_x, const std::uint64_t& p_l, const std::uint64_t& p_r)
    {
        p_h *= p_x;
        p_h ^= (p_h << p_l) | (p_h >> (64 - p_l));
        p_h ^= (p_h >> p_r) | (p_h << (64 - p_r));
    }

    void hash_block_d(std::uint64_t& p_h, const std::vector<std::uint64_t>& p_coeffs)
    {
        for (size_t i = 0; i + 2 < p_coeffs.size(); i += 3)
        {
            hash_block_d(p_h, p_coeffs[i], p_coeffs[i+1], p_coeffs[i+2]);
        }
    }

    class xlr_hash
    {
    public:
        static std::uint64_t hash(const std::uint64_t& p_seed, const std::uint64_t& p_x)
        {
            //3b7d530f5d27b987
            //317371dc620c3879
            //1fce11d96fa69a03
            //fc3c906111de6791
            std::uint64_t h = p_seed ^ p_x;
            hash_block<0x3b7d530f5d27b987ULL,15,13,0xfc3c906111de6791ULL,20,19>(h);
            return h;
        }

        xlr_hash(const std::vector<std::uint64_t>& p_coeffs)
            : m_coeffs(p_coeffs)
        {
        }

        std::uint64_t operator()(const std::uint64_t& p_seed, const std::uint64_t& p_x) const
        {
            std::uint64_t h = p_seed ^ p_x;
            hash_block_d(h, m_coeffs);
            return h;
        }

        std::uint64_t operator()(const std::uint64_t& p_seed, const std::string& p_x) const
        {
            std::uint64_t h = p_seed;
            size_t j = 0;
            std::uint64_t y = 0;
            for (size_t i = 0; i < p_x.size(); ++i)
            {
                y = (y << 8) | p_x[i];
                if (++j == 8)
                {
                    h ^= y;
                    hash_block_d(h, m_coeffs);
                    y = 0;
                    j = 0;
                }
            }
            if (j != 0)
            {
                h ^= y;
                hash_block_d(h, m_coeffs);
            }
            return h;
        }

    private:
        const std::vector<std::uint64_t> m_coeffs;
    };
}
// namespace varoom

#endif // VAROOM_XLR_HASH_HPP

