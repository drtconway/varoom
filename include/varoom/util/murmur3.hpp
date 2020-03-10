#ifndef VAROOM_UTIL_MURMER3_HPP
#define VAROOM_UTIL_MURMER3_HPP

namespace varoom
{
    // Implementation based on Wikipedia page.
    //
    class murmur3
    {
    private:
        static constexpr uint32_t c1 = 0xcc9e2d51U;
        static constexpr uint32_t c2 = 0x1b873593U;
        static constexpr uint32_t r1 = 15;
        static constexpr uint32_t r2 = 13;
        static constexpr uint32_t m = 5;
        static constexpr uint32_t n = 0xe6546b64U;

    public:
        murmur3(uint32_t p_seed)
            : h(p_seed), l(0)
        {
        }

        murmur3& update(uint64_t x)
        {
            uint32_t hi = static_cast<uint32_t>(x >> 32);
            uint32_t lo = static_cast<uint32_t>(x);
            return update(hi).update(lo);
        }

        murmur3& update(uint32_t k)
        {
            k *= c1;
            k = rotl(k, r1);
            k *= c2;

            h ^= k;
            h = rotl(k, r2);
            h = (h * m) + n;

            l += 4;

            return *this;
        }

        uint32_t operator()() const
        {
            uint32_t r = h ^ l;
            r ^= (r >> 16);
            r *= 0x85ebca6bU;
            r ^= (r >> 13);
            r *= 0xc2b2ae35U;
            r ^= (r >> 16);
            return r;
        }

    private:
        static inline uint32_t rotl(uint32_t x, uint32_t r)
        {
            return (x << r) ^ (x >> (32 - r));
        }

        static inline uint32_t rotr(uint32_t x, uint32_t r)
        {
            return (x >> r) ^ (x << (32 - r));
        }

        uint32_t h;
        uint64_t l;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_MURMER3_HPP
