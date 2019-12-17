#ifndef VAROOM_UTIL_GDBM_HPP
#define VAROOM_UTIL_GDBM_HPP

#include <string>
#include <stdexcept>
#include <cstdlib>
#include <gdbm.h>

namespace varoom
{
    class gdbm
    {
    public:
        using datum = std::pair<char*,size_t>;

        template <typename T>
        static void make_datum(const T& p_itm, gdbm::datum& p_datum)
        {
            p_datum.first = reinterpret_cast<const char*>(&p_itm);
            p_datum.second = sizeof(T);
        }

        template <typename T>
        static gdbm::datum make_datum(const T& p_itm)
        {
            const char* ptr = reinterpret_cast<const char*>(&p_itm);
            gdbm::datum d;
            d.first = const_cast<char*>(ptr);
            d.second = sizeof(T);
            return d;
        }

        static void str_datum(const std::string& p_itm, gdbm::datum& p_datum)
        {
            p_datum.first = const_cast<char*>(p_itm.c_str());
            p_datum.second = p_itm.size();
        }

        static gdbm::datum str_datum(const std::string& p_itm)
        {
            gdbm::datum d;
            d.first = const_cast<char*>(p_itm.c_str());
            d.second = p_itm.size();
            return d;
        }

        template <typename T>
        static T datum_cast(const gdbm::datum& p_datum)
        {
            if (p_datum.second != sizeof(T))
            {
                throw std::runtime_error("cannot perform memory_cast");
            }
            const T* res = reinterpret_cast<const T*>(p_datum.first);
            return *res;
        }

        template <typename T>
        static void datum_cast(const gdbm::datum& p_datum, T& p_res)
        {
            if (p_datum.second != sizeof(T))
            {
                throw std::runtime_error("cannot perform memory_cast");
            }
            p_res = *reinterpret_cast<const T*>(p_datum.first);
        }

        struct owning_datum : datum
        {
            owning_datum()
            {
                first = NULL;
                second = 0;
            }

            owning_datum(char* p_data, size_t p_size)
            {
                first = p_data;
                second = p_size;
            }

            owning_datum(const owning_datum& p_other) = delete;
            owning_datum& operator=(const owning_datum& p_other) = delete;

            owning_datum(owning_datum&& p_other)
            {
                first = p_other.first;
                second = p_other.second;
                p_other.first = NULL;
                p_other.second = 0;
            }

            owning_datum& operator=(owning_datum&& p_other)
            {
                if (first)
                {
                    free(first);
                }
                first = p_other.first;
                second = p_other.second;
                p_other.first = NULL;
                p_other.second = 0;
                return *this;
            }

            std::string str() const
            {
                return std::string(first, second);
            }

            template <typename T>
            const T& as() const
            {
                if (second != sizeof(T))
                {
                    throw std::runtime_error("as: type size mismatch");
                }
                return *reinterpret_cast<const T*>(first);
            }

            ~owning_datum()
            {
                if (first)
                {
                    free(first);
                }
            }
        };

        enum flags { read_only, read_write, trunc };
        gdbm(const std::string& p_name, flags p_flags)
            : m_file(NULL)
        {
            int f = 0;
            switch (p_flags)
            {
                case read_only:
                {
                    f = GDBM_READER;
                    break;
                }
                case read_write:
                {
                    f = GDBM_WRCREAT;
                    break;
                }
                case trunc:
                {
                    f = GDBM_NEWDB;
                    break;
                }
            }

            m_file = gdbm_open(p_name.c_str(), 4096, f, 0666, 0);
            if (!m_file)
            {
                throw std::runtime_error("unable to create gdbm object");
            }
        }

        size_t size()
        {
            gdbm_count_t cx = 0;
            int err = gdbm_count(m_file, &cx);
            if (err < 0)
            {
                throw std::runtime_error("unable to determine gdbm size");
            }
            return cx;
        }

        size_t count(const std::string& p_key)
        {
            ::datum key = { const_cast<char*>(p_key.c_str()), static_cast<int>(p_key.size()) };
            return gdbm_exists(m_file, key);
        }
        
        bool get(const gdbm::datum& p_key, gdbm::owning_datum& p_val)
        {
            ::datum key;
            to_raw_datum(p_key, key);
            ::datum val = gdbm_fetch(m_file, key);
            if (val.dptr == NULL)
            {
                return false;
            }
            else
            {
                from_raw_datum(val, p_val);
                return true;
            }
        }

        bool get(const std::string& p_key, std::string& p_val)
        {
            gdbm::datum key;
            str_datum(p_key, key);
            gdbm::owning_datum val;
            if (get(key, val))
            {
                p_val = std::string(val.first, val.second);
                return true;
            }
            else
            {
                return false;
            }
        }

        gdbm::owning_datum get(const std::string& p_key)
        {
            gdbm::datum key;
            str_datum(p_key, key);
            gdbm::owning_datum val;
            if (get(key, val))
            {
                return val;
            }
            else
            {
                throw std::runtime_error("no such key");
            }
        }

        void put(const std::string& p_key, const std::string& p_val)
        {
            gdbm::datum key;
            str_datum(p_key, key);
            gdbm::datum val;
            str_datum(p_val, val);
            put(key, val);
        }

        void put(const gdbm::datum& p_key, const gdbm::datum& p_val)
        {
            ::datum key;
            to_raw_datum(p_key, key);
            ::datum val;
            to_raw_datum(p_val, val);
            gdbm_store(m_file, key, val, GDBM_REPLACE);
        }

        gdbm::owning_datum first()
        {
            ::datum key = gdbm_firstkey(m_file);
            return gdbm::owning_datum(key.dptr, key.dsize);
        }

        bool next(gdbm::owning_datum& p_key)
        {
            if (p_key.first == NULL)
            {
                ::datum key = gdbm_firstkey(m_file);
                if (key.dptr == NULL)
                {
                    return false;
                }
                p_key = gdbm::owning_datum(key.dptr, key.dsize);
                return true;
            }

            ::datum key { p_key.first, static_cast<int>(p_key.second) };
            ::datum nxt = gdbm_nextkey(m_file, key);
            if (nxt.dptr == NULL)
            {
                return false;
            }
            p_key = gdbm::owning_datum(nxt.dptr, nxt.dsize);
            return true;
        }

        ~gdbm()
        {
            if (m_file)
            {
                gdbm_close(m_file);
                m_file = NULL;
            }
        }

    private:
        static void to_raw_datum(const gdbm::datum& p_item, ::datum& p_raw)
        {
            p_raw.dptr = const_cast<char*>(p_item.first);
            p_raw.dsize = static_cast<int>(p_item.second);
        }

        static void from_raw_datum(::datum& p_raw, gdbm::owning_datum& p_item)
        {
            p_item.first = p_raw.dptr;
            p_item.second = p_raw.dsize;
        }

        GDBM_FILE m_file;
    };

}
// namespace varoom

#endif // VAROOM_UTIL_GDBM_HPP
