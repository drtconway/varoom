#ifndef VAROOM_UTIL_TSV_HPP
#define VAROOM_UTIL_TSV_HPP

#ifndef VAROOM_UTIL_SUBTEXT_HPP
#include "varoom/util/subtext.hpp"
#endif

#include <istream>
#include <ostream>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <nlohmann/json.hpp>

namespace varoom
{
    typedef nlohmann::json json;

    class tsv_row
    {
    public:
        tsv_row(const std::unordered_map<std::string,size_t>& p_hdr,
                const std::vector<subtext>& p_row)
            : m_hdr(p_hdr), m_row(p_row)
        {
        }

        size_t size() const
        {
            return m_row.size();
        }

        const subtext& operator[](const size_t& p_idx) const
        {
            return m_row[p_idx];
        }

        const subtext& operator[](const char*& p_fld) const
        {
            auto itr = m_hdr.find(p_fld);
            if (itr == m_hdr.end())
            {
                throw std::runtime_error("field name not in header");
            }
            size_t idx = itr->second;
            return (*this)[idx];
        }

        const subtext& operator[](const std::string& p_fld) const
        {
            auto itr = m_hdr.find(p_fld);
            if (itr == m_hdr.end())
            {
                throw std::runtime_error("field name not in header");
            }
            size_t idx = itr->second;
            return (*this)[idx];
        }

        template <typename T>
        T get(const size_t& p_idx) const
        {
            const subtext& s = (*this)[p_idx];
            return boost::lexical_cast<T>(boost::make_iterator_range(s.first, s.second));
        }

        template <typename T>
        T get(const char*& p_fld) const
        {
            const subtext& s = (*this)[p_fld];
            return boost::lexical_cast<T>(boost::make_iterator_range(s.first, s.second));
        }

        template <typename T>
        T get(const std::string& p_fld) const
        {
            const subtext& s = (*this)[p_fld];
            return boost::lexical_cast<T>(boost::make_iterator_range(s.first, s.second));
        }

    private:
        const std::unordered_map<std::string,size_t>& m_hdr;
        const std::vector<subtext>& m_row;
    };

    class tsv_reader
    {
    public:
        tsv_reader(std::istream& p_in, bool p_header = true)
            : m_in(p_in), m_curr(m_hdr, m_parts), m_more(true)
        {
            next();
            if (p_header && more())
            {
                make_header();
                next();
            }
        }

        const std::unordered_map<std::string,size_t>& header() const
        {
            return m_hdr;
        }
    
        bool more() const
        {
            return m_more;
        }

        const tsv_row& operator*() const
        {
            return m_curr;
        }

        void operator++()
        {
            next();
        }

    private:
        void make_header()
        {
            for (size_t i = 0; i < m_parts.size(); ++i)
            {
                m_hdr[std::string(m_parts[i].first, m_parts[i].second)] = i;
            }
        }

        void next()
        {
            while (m_more)
            {
                if (!std::getline(m_in, m_line))
                {
                    m_more = false;
                    return;
                }
                if (starts_with(m_line, '#'))
                {
                    // comments
                    continue;
                }
                break;
            }

            subtext st(m_line);
            st.split('\t', m_parts);
        }

        static bool starts_with(const std::string& p_str, char p_ch)
        {
            return p_str.size() > 0 && p_str.front() == p_ch;
        }

        std::istream& m_in;
        std::string m_line;
        std::unordered_map<std::string,size_t> m_hdr;
        std::vector<subtext> m_parts;
        tsv_row m_curr;
        bool m_more;
    };

    using tsv_formatter = std::string(*)(const json& p_value);

    class tsv_writer
    {
    public:
        tsv_writer(std::ostream& p_out)
            : m_out(p_out)
        {
        }

        void add(const json& p_row)
        {
            if (m_formatters.size())
            {
                size_t i = 0;
                for (auto itr = p_row.begin(); itr != p_row.end(); ++itr, ++i)
                {
                    if (i > 0)
                    {
                        m_out << '\t';
                    }
                    m_out << m_formatters[i](*itr);
                }
                m_out << std::endl;
            }
            else
            {
                size_t i = 0;
                for (auto itr = p_row.begin(); itr != p_row.end(); ++itr, ++i)
                {
                    if (i > 0)
                    {
                        m_out << '\t';
                    }
                    m_out << (*itr);
                }
                m_out << std::endl;
            }
        }

        void end()
        {
        }

    private:
        std::ostream& m_out;
        std::vector<tsv_formatter> m_formatters;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_TSV_HPP
