#ifndef VAROOM_UTIL_SUBTEXT_HPP
#define VAROOM_UTIL_SUBTEXT_HPP

#include <string>
#include <utility>
#include <vector>

namespace varoom
{
    class subtext : public std::pair<std::string::const_iterator,std::string::const_iterator>
    {
    public:
        subtext(const std::string& p_str)
        {
            first = p_str.begin();
            second = p_str.end();
        }

        subtext(const std::string::const_iterator& p_begin, const std::string::const_iterator& p_end)
        {
            first = p_begin;
            second = p_end;
        }

        size_t size() const
        {
            return second - first;
        }

        bool operator==(const std::string& p_other) const
        {
            return *this == subtext(p_other.begin(), p_other.end());
        }

        bool operator==(const subtext& p_other) const
        {
            std::string::const_iterator i = first;
            std::string::const_iterator j = p_other.first;
            while (i != second && j != p_other.second)
            {
                if (*i++ != *j++)
                {
                    return false;
                }
            }
            return i == second && j == p_other.second;
        }

        operator std::string() const
        {
            return std::string(first, second);
        }

        void split(char p_x, std::vector<subtext>& p_res, bool p_include_empty = true) const
        {
            p_res.clear();
            std::string::const_iterator prev = first;
            while (true)
            {
                std::string::const_iterator next = find(prev, p_x);
                if (p_include_empty || next != prev)
                {
                    p_res.push_back(subtext(prev, next));
                }
                if (next == second)
                {
                    break;
                }
                prev = next + 1;
            }
        }

    private:

        std::string::const_iterator find(char p_x) const
        {
            return find(first, p_x);
        }

        std::string::const_iterator find(const std::string::const_iterator& p_from, char p_x) const
        {
            for (auto itr = p_from; itr != second; ++itr)
            {
                if (*itr == p_x)
                {
                    return itr;
                }
            }
            return second;
        }
    };
}
// namespace varoom

#endif // VAROOM_UTIL_SUBTEXT_HPP
