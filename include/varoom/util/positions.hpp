#ifndef VAROOM_UTIL_POSITIONS_HPP
#define VAROOM_UTIL_POSITIONS_HPP

#ifndef VAROOM_UTIL_STRONG_TYPEDEF_HPP
#include "varoom/util/strong_typedef.hpp"
#endif

#include <string>
#include <boost/flyweight.hpp>

namespace varoom
{
    typedef boost::flyweight<std::string> accession;

    struct position : varoom::strong_typedef<position, uint64_t>, varoom::integer_arithmetic<position>
    {
        using strong_typedef::strong_typedef;
    };

    struct location
    {
        accession chr;
        position pos;

        location()
            : chr(""), pos(0)
        {
        }

        location(const std::string& p_str)
        {
            size_t i = 0;
            while (i < p_str.size() && p_str[i] != ':')
            {
                pos.push_back(p_str[i]);
                ++i;
            }
            if (i == p_str.size() || p_str[i] != ':')
            {
                throw std::runtime_error("unable to parse location - expected ':'");
            }
            ++i
            if (i == p_str.size())
            {
                throw std::runtime_error("unable to parse location - expected digits");
            }
            uint64_t n = 0;
            while (i < p_str.size())
            {
                char c = p_str[i];
                if ('0' <= c && c <= '9')
                {
                    n = 10*n + (c - '0');
                }
                else
                {
                    throw std::runtime_error("unable to parse location - expected digit");
                }
                ++i;
            }
            pos = position(n);
        }
    };

    struct locus
    {
        accession chr;
        position start;
        position end;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_POSITIONS_HPP
