#ifndef VAROOM_UTIL_CASKET_HPP
#define VAROOM_UTIL_CASKET_HPP

#ifndef VAROOM_UTIL_FILES_HPP
#include "varoom/util/files.hpp"
#endif

#include <vector>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/restrict.hpp>
#include <boost/iostreams/stream.hpp>

namespace varoom
{
    namespace detail
    {
        class bytes_input_file_holder : public input_file_holder
        {
        public:
            bytes_input_file_holder(const std::string& p_src)
                : m_device(&p_src[0], p_src.size()), m_stream(m_device)
            {
            }

            std::istream& operator*()
            {
                return m_stream;
            }

        private:
            boost::iostreams::array_source m_device;
            boost::iostreams::stream<boost::iostreams::array_source> m_stream;
        };

        class caching_bytes_input_file_holder : public input_file_holder
        {
        public:
            caching_bytes_input_file_holder(std::string&& p_src)
                : m_src(p_src), m_device(&m_src[0], m_src.size()), m_stream(m_device)
            {
            }

            std::istream& operator*()
            {
                return m_stream;
            }

        private:
            std::string m_src;
            boost::iostreams::array_source m_device;
            boost::iostreams::stream<boost::iostreams::array_source> m_stream;
        };

        class bytes_output_file_holder : public output_file_holder
        {
        public:
            bytes_output_file_holder(std::function<void(const std::vector<char>&)> p_save)
                : m_save(p_save), m_device(m_bytes), m_stream(m_device)
            {
            }

            ~bytes_output_file_holder()
            {
                (**this).flush();
                m_save(m_bytes);
            }

            std::ostream& operator*()
            {
                return m_stream;
            }

        private:
            std::function<void(const std::vector<char>&)> m_save;
            std::vector<char> m_bytes;
            boost::iostreams::back_insert_device<std::vector<char>> m_device;
            boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char>>> m_stream;
        };

        class ranged_input_file_holder : public input_file_holder
        {
        public:
            ranged_input_file_holder(std::istream& p_in, size_t p_off, size_t p_len)
                : m_stream(boost::iostreams::restrict(p_in, p_off, p_len))
            {
            }

            std::istream& operator*()
            {
                return m_stream;
            }

        private:
            boost::iostreams::filtering_istream m_stream;
        };

    }

    class casket
    {
    public:
        using offset_and_length = std::pair<uint64_t,uint64_t>;
        using toc_type = std::map<std::string,offset_and_length>;

        casket(std::istream& p_in)
            : m_src(&p_in), m_dst(NULL), m_dst_pos(0)
        {
            m_src->seekg(0, std::ios::end);
            m_src_len = m_src->tellg();
            read_toc();
        }

        casket(std::ostream& p_out)
            : m_src(NULL), m_src_len(0), m_dst(&p_out), m_dst_pos(0)
        {
        }

        ~casket()
        {
            if (m_dst)
            {
                write_toc();
            }
        }

        std::vector<std::string> contents() const
        {
            std::vector<std::string> res;
            contents(res);
            return res;
        }

        void contents(std::vector<std::string>& p_res) const
        {
            p_res.clear();
            for (auto itr = m_toc.begin(); itr != m_toc.end(); ++itr)
            {
                p_res.push_back(itr->first);
            }
        }

        input_file_holder_ptr in(const std::string& p_name)
        {
            require_src();

            auto itr = m_toc.find(p_name);
            if (itr == m_toc.end())
            {
                throw std::runtime_error("no such file in casket");
            }
            uint64_t p = itr->second.first;
            uint64_t l = itr->second.second;
            m_src->seekg(p, std::ios_base::beg);
            std::string s;
            s.resize(l);
            m_src->read(&s[0], l);
            return input_file_holder_ptr(new detail::caching_bytes_input_file_holder(std::move(s)));
        }

        output_file_holder_ptr out(const std::string& p_name)
        {
            std::function<void(const std::vector<char>&)> f = [this,p_name](const std::vector<char>& p_bytes) mutable {
                write_file(p_name, p_bytes);
            };
            return output_file_holder_ptr(new detail::bytes_output_file_holder(f));
        }

        void with(const std::string& p_name, std::function<void(std::istream&)> p_func)
        {
            require_src();

            auto itr = m_toc.find(p_name);
            if (itr == m_toc.end())
            {
                throw std::runtime_error("no such file in casket");
            }
            uint64_t p = itr->second.first;
            uint64_t l = itr->second.second;
            m_src->seekg(0, std::ios_base::beg);
            detail::ranged_input_file_holder i(*m_src, p, l);
            p_func(*i);
        }

        void with(const std::string& p_name, std::function<void(std::ostream&)> p_func)
        {
            require_src();

            auto itr = m_toc.find(p_name);
            if (itr == m_toc.end())
            {
                throw std::runtime_error("no such file in casket");
            }
            std::vector<char> b;
            boost::iostreams::back_insert_device<std::vector<char>> d(b);
            boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char>>> s(d);
            p_func(s);
            s.flush();
            write_file(p_name, b);
        }

    private:
        void read_toc()
        {
            require_src();

            // Seek to the end and read the offset of the TOC.
            m_src->seekg(m_src_len - sizeof(uint64_t));
            uint64_t z = read_pod<uint64_t>();

            m_src->seekg(z, std::ios::beg);
            // Get the number of entries.
            uint64_t n = read_pod<uint64_t>();
            for (size_t i = 0; i < n; ++i)
            {
                // Position & Length
                offset_and_length v;
                v.first = read_pod<uint64_t>();
                v.second = read_pod<uint64_t>();
                std::string nm = read_str();
                m_toc[nm] = v;
            }
        }

        void write_toc()
        {
            require_dst();

            uint64_t toc_pos = m_dst_pos;

            uint64_t toc_len = m_toc.size();
            write_pod(toc_len);
            for (auto itr = m_toc.begin(); itr != m_toc.end(); ++itr)
            {
                write_pod(itr->second.first);
                write_pod(itr->second.second);
                write_str(itr->first);
            }
            write_pod(toc_pos);
        }

        void write_file(const std::string& p_name, const std::vector<char>& p_bytes)
        {
            offset_and_length v;
            v.first = m_dst_pos;
            v.second = p_bytes.size();
            m_toc[p_name] = v;

            m_dst->write(p_bytes.data(), p_bytes.size());
            m_dst_pos += p_bytes.size();
        }

        template <typename T>
        T read_pod()
        {
            static_assert(std::is_pod<T>::value);
            require_src();
            T x;
            m_src->read(reinterpret_cast<char*>(&x), sizeof(x));
            return x;
        }

        template <typename T>
        void write_pod(const T& p_x)
        {
            static_assert(std::is_pod<T>::value);
            require_dst();
            m_dst->write(reinterpret_cast<const char*>(&p_x), sizeof(T));
            m_dst_pos += sizeof(T);
        }

        void read_str(std::string& p_str)
        {
            require_src();
            uint64_t n = read_pod<uint64_t>();
            p_str.resize(n);
            m_src->read(&p_str[0], n);
        }

        std::string read_str()
        {
            std::string r;
            read_str(r);
            return r;
        }

        void write_str(const std::string& p_str)
        {
            require_dst();
            uint64_t l = p_str.size();
            write_pod(l);
            m_dst->write(&p_str[0], l);
            m_dst_pos += l;
        }

        void require_src()
        {
            if (!m_src)
            {
                throw std::runtime_error("cannot invoke method without input stream");
            }
        }

        void require_dst()
        {
            if (!m_dst)
            {
                throw std::runtime_error("cannot invoke method without output stream");
            }
        }

        std::istream* m_src;
        size_t m_src_len;
        std::ostream* m_dst;
        size_t m_dst_pos;
        toc_type m_toc;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_CASKET_HPP
