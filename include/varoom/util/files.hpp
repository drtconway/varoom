#ifndef VAROOM_UTIL_FILES_HPP
#define VAROOM_UTIL_FILES_HPP

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace varoom
{
    class input_file_holder
    {
    public:
        virtual std::istream& operator*() = 0;

        virtual ~input_file_holder() {}
    };
    typedef std::shared_ptr<input_file_holder> input_file_holder_ptr;

    namespace detail
    {
        class cin_file_holder : public input_file_holder
        {
        public:
            std::istream& operator*()
            {
                return std::cin;
            };
        };

        class plain_file_holder : public input_file_holder
        {
        public:
            plain_file_holder(const std::string& p_name)
                : m_file(p_name)
            {
            }

            std::istream& operator*()
            {
                return m_file;
            }

        private:
            std::ifstream m_file;
        };

        class gzip_file_holder : public input_file_holder
        {
        public:
            gzip_file_holder(const std::string& p_name)
                : m_file(p_name, std::ios_base::in | std::ios_base::binary)
            {
                m_filter.push(m_gzip);
                m_filter.push(m_file);
            }

            std::istream& operator*()
            {
                return m_filter;
            }

        private:
            std::ifstream m_file;
            boost::iostreams::gzip_decompressor m_gzip;
            boost::iostreams::filtering_istream m_filter;
        };
    }
    // namespace detail

    class files
    {
    public:
        static input_file_holder_ptr in(const std::string& p_name)
        {
            if (p_name == "-")
            {
                return input_file_holder_ptr(new detail::cin_file_holder);
            }
            if (ends_with(p_name, ".gz"))
            {
                return input_file_holder_ptr(new detail::gzip_file_holder(p_name));
            }
            return input_file_holder_ptr(new detail::plain_file_holder(p_name));
        }

    private:
        static bool ends_with(const std::string& p_str, const std::string& p_suffix)
        {
            auto s = p_str.rbegin();
            auto t = p_suffix.rend();
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

#endif // VAROOM_UTIL_FILES_HPP
