#ifndef VAROOM_UTIL_TEXT_HPP
#define VAROOM_UTIL_TEXT_HPP

#include <initializer_list>

namespace varoom
{
    class text
    {
    public:
        static char to_upper(char p_ch)
        {
            if ('a' <= p_ch && p_ch <= 'z')
            {
                return 'A' + (p_ch - 'a');
            }
            return p_ch;
        }

        static std::string tabs(std::initializer_list<const char*> p_parts)
        {
            std::string s;
            for (auto i = p_parts.begin(); i != p_parts.end(); ++i)
            {
                if (s.size() > 0)
                {
                    s.push_back('\t');
                }
                s.insert(s.size(), *i);
            }
            return s;
        }

        static bool starts_with(const std::string& p_str, char p_ch)
        {
            return p_str.size() > 0 && p_str.front() == p_ch;
        }

        static bool starts_with(const std::string& p_str, const std::string& p_prefix)
        {
            auto s = p_str.begin();
            auto t = p_prefix.begin();
            while (s != p_str.end() && t != p_prefix.end())
            {
                if (*s != *t)
                {
                    return false;
                }
                ++s;
                ++t;
            }
            return t == p_prefix.end();
        }

        static bool ends_with(const std::string& p_str, const std::string& p_suffix)
        {
            auto s = p_str.rbegin();
            auto t = p_suffix.rbegin();
            while (s != p_str.rend() && t != p_suffix.rend())
            {
                if (*s != *t)
                {
                    return false;
                }
                ++s;
                ++t;
            }
            return t == p_suffix.rend();
        }
    };
}
// namespace varoom

#endif // VAROOM_UTIL_TEXT_HPP
