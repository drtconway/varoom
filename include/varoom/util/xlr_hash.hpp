#ifndef VAROOM_XLR_HASH_HPP
#define VAROOM_XLR_HASH_HPP

#include <cstdint>

namespace varoom
{
    template<unsigned long long X, unsigned long long L, unsigned long long R>
    void hash_block(std::uint64_t& p_h)
    {
        p_h ^= X;
        p_h ^= (p_h << L) | (p_h >> R);
    }

    template<unsigned long long X, unsigned long long L, unsigned long long R, unsigned long long Y, unsigned long long... More>
    void hash_block(std::uint64_t& p_h)
    {
        hash_block<X,L,R>(p_h);
        hash_block<Y,More...>(p_h);
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
    };
}
// namespace varoom

#endif // VAROOM_XLR_HASH_HPP

