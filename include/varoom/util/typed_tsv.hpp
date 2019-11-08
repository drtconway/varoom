#ifndef VAROOM_UTIL_TYPED_TSV_HPP
#define VAROOM_UTIL_TYPED_TSV_HPP

#ifndef VAROOM_UTIL_TSV_HPP
#include "varoom/util/tsv.hpp"
#endif

#include <boost/any.hpp>

namespace varoom
{
    typedef boost::any tsv_column_value;
    typedef std::vector<tsv_column_value> typed_tsv_row;

    class tsv_column_type;
    typedef std::shared_ptr<tsv_column_type> tsv_column_type_ptr;

    class tsv_column_type
    {
    public:
        tsv_column_type(const std::string& p_name)
            : m_name(p_name)
        {
        }

        virtual ~tsv_column_type() {}

        const std::string& name() const
        {
            return m_name;
        }

        virtual tsv_column_value make(const subtext& p_txt) const = 0;

        virtual void unmake(const tsv_column_value& p_val, std::string& p_str) const = 0;

        static bool add(tsv_column_type_ptr p_type)
        {
            types()[p_type->name()] = p_type;
            return true;
        }

        static const tsv_column_type_ptr& get(const std::string& p_name)
        {
            auto itr = types().find(p_name);
            if (itr == types().end())
            {
                throw std::runtime_error("no such registered type");
            }
            return itr->second;
        }

    private:
        static std::unordered_map<std::string,tsv_column_type_ptr>& types()
        {
            static std::unordered_map<std::string,tsv_column_type_ptr> m;
            return m;
        }

        const std::string m_name;
    };

    class typed_tsv_reader
    {
    public:
        typed_tsv_reader(std::istream& p_in, const std::vector<std::string>& p_types, bool p_header = true)
            : m_tsv(p_in, p_header), m_type_names(p_types)
        {
            for (size_t i = 0; i < m_type_names.size(); ++i)
            {
                m_types.push_back(tsv_column_type::get(m_type_names[i]));
            }
            if (more())
            {
                load_row();
            }
        }

        bool more() const
        {
            return m_tsv.more();
        }

        const typed_tsv_row& operator*() const
        {
            return m_curr;
        }

        void operator++()
        {
            ++m_tsv;
            load_row();
        }
    private:
        void load_row()
        {
            const tsv_row& row = *m_tsv;
            if (row.size() != m_types.size())
            {
                throw std::runtime_error("mismatch between number of elements in row and number of types");
            }
            m_curr.resize(row.size());
            for (size_t i = 0; i < m_types.size(); ++i)
            {
                m_curr[i] = m_types[i]->make(row[i]);
            }
        }

        tsv_reader m_tsv;
        typed_tsv_row m_curr;
        std::vector<std::string> m_type_names;
        std::vector<tsv_column_type_ptr> m_types;
    };

    namespace detail
    {
        template <typename T>
        class tsv_column_type_atom : public tsv_column_type
        {
        public:
            tsv_column_type_atom(const std::string& p_name)
                : tsv_column_type(p_name)
            {
            }

            virtual tsv_column_value make(const subtext& p_txt) const
            {
                T x = boost::lexical_cast<T>(boost::make_iterator_range(p_txt.first, p_txt.second));
                return tsv_column_value(x);
            }

            virtual void unmake(const tsv_column_value& p_val, std::string& p_str) const
            {
                T x = boost::any_cast<T>(p_val);
                p_str = boost::lexical_cast<std::string>(x);
            }
        };
        bool default_tsv_column_atomic_types[] = {
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_atom<int64_t>("int"))),
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_atom<uint64_t>("uint"))),
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_atom<double>("flt"))),
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_atom<std::string>("str")))
        };

        template <typename T>
        class tsv_column_type_key_value : public tsv_column_type
        {
        public:
            typedef std::pair<std::string,T> key_value_pair;

            tsv_column_type_key_value(const std::string& p_name, const std::string& p_type, char p_sep = '=')
                : tsv_column_type(p_name),
                  m_elem_type(*tsv_column_type::get(p_type)),
                  m_sep(p_sep)
            {
            }

            virtual tsv_column_value make(const subtext& p_txt) const
            {
                std::vector<subtext> parts;
                p_txt.split(m_sep, parts);
                if (parts.size() < 2)
                {
                    throw std::runtime_error("malformed key value pair");
                }
                key_value_pair x;

                x.first = parts[0];

                subtext val_txt(parts[1].first, parts.back().second);
                tsv_column_value val = m_elem_type.make(val_txt);
                x.second = boost::any_cast<T>(val);
                
                return tsv_column_value(x);
            }

            virtual void unmake(const tsv_column_value& p_val, std::string& p_str) const
            {
                const key_value_pair& x = boost::any_cast<const key_value_pair&>(p_val);
                tsv_column_value v(x.second);

                std::string s;
                m_elem_type.unmake(v, s);

                p_str = x.first;
                p_str.push_back(m_sep);
                p_str.insert(p_str.end(), s.begin(), s.end());
            }

        private:
            const tsv_column_type& m_elem_type;
            const char m_sep;
        };

        bool default_tsv_column_key_value_types[] = {
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_key_value<int64_t>("str->int", "int"))),
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_key_value<uint64_t>("str->uint", "uint"))),
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_key_value<double>("str->flt", "flt"))),
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_key_value<std::string>("str->str", "str")))
        };

        template <typename T>
        class tsv_column_type_list : public tsv_column_type
        {
        public:
            typedef std::vector<T> list_type;

            tsv_column_type_list(const std::string& p_name, const std::string& p_type, char p_sep)
                : tsv_column_type(p_name),
                  m_elem_type(*tsv_column_type::get(p_type)),
                  m_sep(p_sep)
            {
            }

            virtual tsv_column_value make(const subtext& p_txt) const
            {
                list_type xs;
                if (p_txt.size() == 0)
                {
                    return tsv_column_value(xs);
                }

                std::vector<subtext> parts;
                p_txt.split(m_sep, parts);
                for (size_t i = 0; i < parts.size(); ++i)
                {
                    tsv_column_value v = m_elem_type.make(parts[i]);
                    xs.push_back(boost::any_cast<T>(v));
                }
                return tsv_column_value(xs);
            }

            virtual void unmake(const tsv_column_value& p_val, std::string& p_str) const
            {
                const list_type& xs = boost::any_cast<const list_type&>(p_val);
                p_str.clear();
                for (size_t i = 0; i < xs.size(); ++i)
                {
                    if (i > 0)
                    {
                        p_str.push_back(m_sep);
                    }
                    tsv_column_value v(xs[i]);
                    std::string s;
                    m_elem_type.unmake(v, s);
                    p_str.insert(p_str.end(), s.begin(), s.end());
                }
            }

        private:
            const tsv_column_type& m_elem_type;
            const char m_sep;
        };
        
        bool default_tsv_column_list_types[] = {
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_list<int64_t>("[int]", "int", ';'))),
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_list<uint64_t>("[uint]", "uint", ';'))),
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_list<double>("[flt]", "flt", ';'))),
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_list<std::string>("[str]", "str", ';'))),
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_list<std::pair<std::string,int64_t>>("[str->int]", "str->int", ';'))),
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_list<std::pair<std::string,uint64_t>>("[str->uint]", "str->uint", ';'))),
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_list<std::pair<std::string,double>>("[str->flt]", "str->flt", ';'))),
            tsv_column_type::add(tsv_column_type_ptr(new tsv_column_type_list<std::pair<std::string,std::string>>("[str->str]", "str->str", ';'))),
        };
    }
    // namespace detail
}
// namespace varoom

#endif // VAROOM_UTIL_TYPED_TSV_HPP
