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
        struct datum
        {
            char* data;
            size_t size;

            datum()
                : data(NULL), size(0)
            {
            }

            datum(char* p_data, size_t p_size)
                : data(p_data), size(p_size)
            {
            }

            datum(const datum& p_other) = delete;
            datum& operator=(const datum& p_other) = delete;

            datum(datum&& p_other)
            {
                data = p_other.data;
                size = p_other.size;
                p_other.data = NULL;
                p_other.size = 0;
            }

            datum& operator=(datum&& p_other)
            {
                if (data)
                {
                    free(data);
                }
                data = p_other.data;
                size = p_other.size;
                p_other.data = NULL;
                p_other.size = 0;
                return *this;
            }

            std::string str() const
            {
                return std::string(data, size);
            }

            template <typename T>
            const T& as() const
            {
                if (size != sizeof(T))
                {
                    throw std::runtime_error("as: type size mismatch");
                }
                return *reinterpret_cast<const T*>(data);
            }

            ~datum()
            {
                if (data)
                {
                    free(data);
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
        
        gdbm::datum get(const std::string& p_key)
        {
            ::datum key = { const_cast<char*>(p_key.c_str()), static_cast<int>(p_key.size()) };
            ::datum val = gdbm_fetch(m_file, key);
            return gdbm::datum(val.dptr, val.dsize);
        }

        void put(const std::string& p_key, const std::string& p_val)
        {
            ::datum key = { const_cast<char*>(p_key.c_str()), static_cast<int>(p_key.size()) };
            ::datum val = { const_cast<char*>(p_val.c_str()), static_cast<int>(p_val.size()) };
            gdbm_store(m_file, key, val, GDBM_REPLACE);
        }

        gdbm::datum first()
        {
            ::datum key = gdbm_firstkey(m_file);
            return gdbm::datum(key.dptr, key.dsize);
        }

        bool next(gdbm::datum& p_key)
        {
            if (p_key.data == NULL)
            {
                ::datum key = gdbm_firstkey(m_file);
                if (key.dptr == NULL)
                {
                    return false;
                }
                p_key = gdbm::datum(key.dptr, key.dsize);
                return true;
            }

            ::datum key { p_key.data, static_cast<int>(p_key.size) };
            ::datum nxt = gdbm_nextkey(m_file, key);
            if (nxt.dptr == NULL)
            {
                return false;
            }
            p_key = gdbm::datum(nxt.dptr, nxt.dsize);
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
        GDBM_FILE m_file;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_GDBM_HPP
