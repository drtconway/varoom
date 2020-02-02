#ifndef VAROOM_UTIL_ELIAS_FANO_HPP
#define VAROOM_UTIL_ELIAS_FANO_HPP

#ifndef VAROOM_UTIL_BITSTREAM_HPP
#include "varoom/util/bitstream.hpp"
#endif

#include <iostream>

namespace varoom
{
    struct elias
    {
        static size_t ilog2c(uint64_t x)
        {
            size_t n = 0;
            while (x > 0)
            {
                x >>= 1;
                n += 1;
            }
            return n;
        }

        static uint64_t rev(uint64_t x)
        {
            x = ((x & 0x5555555555555555ULL) << 1) | ((x & 0xAAAAAAAAAAAAAAAAULL) >> 1);
            x = ((x & 0x3333333333333333ULL) << 2) | ((x & 0xCCCCCCCCCCCCCCCCULL) >> 2);
            x = ((x & 0x0F0F0F0F0F0F0F0FULL) << 4) | ((x & 0xF0F0F0F0F0F0F0F0ULL) >> 4);
            x = ((x & 0x00FF00FF00FF00FFULL) << 8) | ((x & 0xFF00FF00FF00FF00ULL) >> 8);
            x = ((x & 0x0000FFFF0000FFFFULL) <<16) | ((x & 0xFFFF0000FFFF0000ULL) >>16);
            x = ((x & 0x00000000FFFFFFFFULL) <<32) | ((x & 0xFFFFFFFF00000000ULL) >>32);
            return x;
        }

        static uint64_t rev(size_t n, uint64_t x)
        {
            x = rev(x);
            return x >> (64 - n);
        }

        static varoom::bitstream& encode(varoom::bitstream& p_bits, uint64_t x)
        {
            //assert(x >= 1);

            size_t n = ilog2c(x);
            p_bits.push_back(n - 1, 0ULL);
            p_bits.push_back(n, rev(n, x));

            return p_bits;
        }

        template <typename Bits>
        static size_t count_leading_zeros(Bits& p_bits)
        {
            size_t n = 0;
            while (*p_bits == 0)
            {
                ++n;
                ++p_bits;
            }
            return n;
        }

        template <typename Bits>
        static uint64_t decode(Bits& p_bits)
        {
            size_t n = count_leading_zeros(p_bits) + 1;
            uint64_t x = rev(n, p_bits.pop_front(n));
            return x;
        }
    };

    struct golomb
    {
        static varoom::bitstream& encode(varoom::bitstream& p_bits, uint64_t x)
        {
            return elias::encode(p_bits, x + 1);
        }

        template <typename Bits>
        static uint64_t decode(Bits& p_bits)
        {
            return elias::decode(p_bits) - 1;
        }
    };
}
// namespace varoom


#endif // VAROOM_UTIL_ELIAS_FANO_HPP
