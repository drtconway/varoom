#ifndef VAROOM_UTIL_CODEC8_HPP
#define VAROOM_UTIL_CODEC8_HPP

#include <cstdint>
#include <vector>

namespace varoom
{
    namespace util
    {
        class codec8
        {
        public:
            static size_t encode(const uint64_t& p_x, std::vector<uint8_t>& p_vec)
            {
                size_t n = 0;
                uint64_t x = p_x;
                uint8_t v = 0;
                uint8_t xs[9];
                do
                {
                    xs[n++] = ((x & 127) << 1) | v;
                    x >>= 7;
                    v = 1;
                }
                while (x > 0);
                size_t i = n;
                while (i > 0)
                {
                    p_vec.push_back(xs[--i]);
                }
                return n;
            }

            template <typename IteratorType>
            static uint64_t decode(IteratorType& p_src)
            {
                uint64_t x = 0;
                while (true)
                {
                    uint8_t y = *p_src;
                    ++p_src;
                    x = (x << 7) | (y >> 1);
                    if ((y & 1) == 0)
                    {
                        break;
                    }
                }
                return x;
            }
        };
    }
    // namespace util
}
// namespace varoom


#endif // VAROOM_UTIL_CODEC8_HPP
