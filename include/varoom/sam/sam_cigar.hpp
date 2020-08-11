#ifndef VAROOM_SAM_SAM_CIGAR_HPP
#define VAROOM_SAM_SAM_CIGAR_HPP

#include <regex>
#include <string>
#include <utility>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/range/iterator_range.hpp>

namespace varoom
{
    typedef std::pair<char,uint32_t> cigar_op;

    void decode_cigar(const std::string& p_cigar, std::vector<cigar_op>& p_ops)
    {
        p_ops.clear();
        uint32_t n = 0;
        for (auto itr = p_cigar.begin(); itr != p_cigar.end(); ++itr)
        {
            char c = *itr;
            if ('0' <= c && c <= '9')
            {
                n = 10*n + (c - '0');
            }
            else
            {
                p_ops.push_back(cigar_op(c, n));
                n = 0;
            }
        }

        if (0)
        {
            const std::regex cig_regex("([0-9]+)([DHIMNPSX=])");

            auto cig_begin = std::sregex_iterator(p_cigar.begin(), p_cigar.end(), cig_regex);
            auto cig_end = std::sregex_iterator();
            for (auto itr = cig_begin; itr != cig_end; ++itr)
            {
                std::smatch m = *itr;
                uint32_t n = boost::lexical_cast<uint32_t>(m[1]);
                char c = *(m[2].first);
                p_ops.push_back(cigar_op(c,n));
            }
        }
    }
}
// namespace varoom

#endif // VAROOM_SAM_SAM_CIGAR_HPP
