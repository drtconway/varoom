#ifndef VAROOM_SEQ_LOCUS_ORDERING_HPP
#define VAROOM_SEQ_LOCUS_ORDERING_HPP

#include <unordered_map>

namespace varoom
{
    struct locus_ordering
    {
        static const std::unordered_map<std::string,int>& ord()
        {
            static std::unordered_map<std::string,int> m = {
                {"1", 1}, {"chr1", 1}, {"2", 2}, {"chr2", 2},
                {"3", 3}, {"chr3", 3}, {"4", 4}, {"chr4", 4},
                {"5", 5}, {"chr5", 5}, {"6", 6}, {"chr6", 6},
                {"7", 7}, {"chr7", 7}, {"8", 8}, {"chr8", 8},
                {"9", 9}, {"chr9", 9}, {"10", 10}, {"chr10", 10},
                {"11", 11}, {"chr11", 11}, {"12", 12}, {"chr12", 12},
                {"13", 13}, {"chr13", 13}, {"14", 14}, {"chr14", 14},
                {"15", 15}, {"chr15", 15}, {"16", 16}, {"chr16", 16},
                {"17", 17}, {"chr17", 17}, {"18", 18}, {"chr18", 18},
                {"19", 19}, {"chr19", 19}, {"20", 20}, {"chr20", 20},
                {"21", 21}, {"chr21", 21}, {"22", 22}, {"chr22", 22},
                {"X", 23}, {"chrX", 23}, {"Y", 24}, {"chrY", 24}
            };
            return m;
        }

        static bool equal(const std::string& p_lhs, const std::string& p_rhs)
        {
            auto e = ord().end();
            auto l = ord().find(p_lhs);
            auto r = ord().find(p_rhs);
            if (l == e && r == e)
            {
                return p_lhs == p_rhs;
            }
            if (l == e)
            {
                return false;
            }
            if (r == e)
            {
                return true;
            }
            return l->second == r->second;
        }

        static bool less(const std::string& p_lhs, const std::string& p_rhs)
        {
            auto e = ord().end();
            auto l = ord().find(p_lhs);
            auto r = ord().find(p_rhs);
            if (l == e && r == e)
            {
                return p_lhs < p_rhs;
            }
            if (l == e)
            {
                return false;
            }
            if (r == e)
            {
                return true;
            }
            return l->second < r->second;
        }

    };
}
// namespace varoom

#endif // VAROOM_SEQ_LOCUS_ORDERING_HPP
