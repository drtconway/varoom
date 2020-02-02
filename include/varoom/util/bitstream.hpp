#ifndef VAROOM_UTIL_BITSTREAM_HPP
#define VAROOM_UTIL_BITSTREAM_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace varoom
{
    class bitstream
    {
    public:
        struct iterator
        {
            const bitstream* bits;
            size_t pos;

            iterator(const bitstream* p_bits, size_t p_pos)
                : bits(p_bits), pos(p_pos)
            {
            }

            bool operator*() const
            {
                return (*bits)[pos];
            }

            uint64_t pop_front(size_t n)
            {
                uint64_t x = bits->at(n, pos);
                pos += n;
                return x;
            }

            iterator& operator++()
            {
                ++pos;
                return *this;
            }

            bool operator==(const iterator& other) const
            {
                return bits == other.bits && pos == other.pos;
            }

            bool operator!=(const iterator& other) const
            {
                return !((*this) == other);
            }
        };

        bitstream()
            : m_word_bit(64)
        {
        }

        size_t size() const
        {
            if (m_words.size() == 0)
            {
                return 0;
            }
            return 64*(m_words.size() - 1) + m_word_bit;
        }

        bool operator[](size_t pos) const
        {
            size_t w = pos / 64;
            size_t b = pos % 64;
            uint64_t m = 1ULL << b;
            return (m_words[w] & m) != 0;
        }

        uint64_t at(size_t n, size_t pos) const
        {
            size_t w = pos / 64;
            size_t b = pos % 64;
            if (b + n <= 64)
            {
                uint64_t m = (1ULL << n) - 1;
                return (m_words[w] >> b) & m;
            }

            size_t l = 64 - b;
            size_t h = n - l;
            uint64_t hm = (1ULL << h) - 1;
            return (m_words[w] >> b) | ((m_words[w + 1] & hm) << l);
        }

        iterator begin() const
        {
            return iterator(this, 0);
        }

        iterator end() const
        {
            return iterator(this, size());
        }

        bitstream& push_back(size_t B, uint64_t p_bits)
        {
            if (m_word_bit == 64)
            {
                make_room();
            }

            uint64_t M = (1ULL << B) - 1;
            uint64_t x = p_bits & M;
            if (m_word_bit + B <= 64)
            {
                m_words.back() |= x << m_word_bit;
                m_word_bit += B;
                return *this;
            }

            size_t b = 64 - m_word_bit;
            size_t c = B - b;
            uint64_t bM = (1ULL << b) - 1;

            uint64_t lo = x & bM;
            m_words.back() |= lo << m_word_bit;
            m_word_bit += b;

            make_room();
    
            uint64_t hi = x >> b;
            m_words.back() |= hi;
            m_word_bit += c;

            return *this;
        }

        template <int B>
        bitstream& push_back(uint64_t p_bits)
        {
            return push_back(B, p_bits);
        }

        std::string str() const
        {
            std::string s;
            for (size_t i = 0; i < size(); ++i)
            {
                s.push_back("01"[(*this)[i]]);
            }
            return s;
        }

    private:
        void make_room()
        {
            //assert(m_word_bit == 64);
            m_words.push_back(0);
            m_word_bit = 0;
        }

        std::vector<uint64_t> m_words;
        size_t m_word_bit;
    };

}
// namespace varoom

#endif // VAROOM_UTIL_BITSTREAM_HPP
