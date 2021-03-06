#ifndef VAROOM_UTIL_TABLE_HPP
#define VAROOM_UTIL_TABLE_HPP

#include <initializer_list>
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/range/iterator_range.hpp>

#ifndef VAROOM_UTIL_SUBTEXT_HPP
#include "varoom/util/subtext.hpp"
#endif

#include <iostream>

namespace varoom
{
    namespace detail
    {
        template <typename... Cols>
        struct table_implementation
        {
            using row = std::tuple<Cols...>;
            virtual ~table_implementation() {}
        };
        template <typename... Cols>
        using table_implementation_ptr = std::shared_ptr<table_implementation<Cols...>>;

        template <typename... Cols>
        struct input_table_implementation : table_implementation<Cols...> {};
        template <typename... Cols>
        using input_table_implementation_ptr = std::shared_ptr<input_table_implementation<Cols...>>;

        template <typename... Cols>
        struct output_table_implementation : table_implementation<Cols...> {};
        template <typename... Cols>
        using output_table_implementation_ptr = std::shared_ptr<output_table_implementation<Cols...>>;

        template <typename... Cols>
        class streamed_input_table : public input_table_implementation<Cols...>
        {
        public:
            using row = typename table_implementation<Cols...>::row;

            streamed_input_table(std::istream& p_in, bool p_header = false)
                : m_in(p_in), m_line_number(0), m_header(p_header)
            {
            }

            bool next(row& p_row)
            {
                while (true)
                {
                    ++m_line_number;
                    if (!std::getline(m_in, m_line))
                    {
                        return false;
                    }
                    if (m_line.size() > 0 && m_line[0] == '#')
                    {
                        continue;
                    }
                    if (m_header)
                    {
                        m_header = false;
                        continue;
                    }
                    subtext st(m_line);
                    st.split('\t', m_parts);
                    make_each(m_parts, p_row);
                    return true;
                }
            }

            size_t line_number() const
            {
                return m_line_number;
            }

            template <typename T>
            void make(const subtext& p_stxt, T& p_res)
            {
                p_res = boost::lexical_cast<T>(boost::make_iterator_range(p_stxt.first, p_stxt.second));
            }

            template<std::size_t I = 0>
            typename std::enable_if<I == sizeof...(Cols), void>::type
            make_each(const std::vector<subtext>&, row&)
            {
            }

            template<std::size_t I = 0>
            typename std::enable_if<I < sizeof...(Cols), void>::type
            make_each(const std::vector<subtext>& p_parts, row& p_row)
            {
                using T_I = typename std::tuple_element<I, row>::type;
                const subtext& s = p_parts[I];
                std::get<I>(p_row) = boost::lexical_cast<T_I>(boost::make_iterator_range(s.first, s.second));
                make_each<I + 1>(p_parts, p_row);
            }

        private:
            std::istream& m_in;
            size_t m_line_number;
            bool m_header;
            std::string m_line;
            std::vector<subtext> m_parts;
        };

        template <typename... Cols>
        class streamed_output_table : public output_table_implementation<Cols...>
        {
        public:
            using row = typename table_implementation<Cols...>::row;

            streamed_output_table(std::ostream& p_out)
                : m_out(p_out)
            {
            }

            streamed_output_table(std::ostream& p_out, std::initializer_list<std::string> p_labels)
                : m_out(p_out)
            {
                //Dang! p_labels.size() didn't become a constexpr till c++14.
                //static_assert(p_labels.size() == std::tuple_size<row>::value);
                size_t n = 0;
                for (auto i = p_labels.begin(); i != p_labels.end(); ++i, ++n)
                {
                    // nothing - we're just counting.
                }
                if (n != std::tuple_size<row>::value)
                {
                    throw std::runtime_error("wrong number of arguments given");
                }
                n = 0;
                for (auto i = p_labels.begin(); i != p_labels.end(); ++i, ++n)
                {
                    if (n > 0)
                    {
                        m_out << '\t';
                    }
                    m_out << *i;
                }
                m_out << std::endl;
            }

            streamed_output_table& operator<<(const row& p_row)
            {
                output_each(p_row);
                return *this;
            }
            
            template<std::size_t I = 0>
            typename std::enable_if<I == sizeof...(Cols), void>::type
            output_each(const row&)
            {
                m_out << std::endl;
            }

            template<std::size_t I = 0>
            typename std::enable_if<I < sizeof...(Cols), void>::type
            output_each(const row& p_row)
            {
                //using T_I = typename std::tuple_element<I, row>::type;
                if (I > 0)
                {
                    m_out << '\t';
                }
                m_out << std::get<I>(p_row);
                output_each<I + 1>(p_row);
            }

        private:
            std::ostream& m_out;
        };

        template <typename... Cols>
        class inmemory_table : public std::vector<std::tuple<Cols...>>, table_implementation<Cols...>
        {
        public:
            using row = typename table_implementation<Cols...>::row;
        };

        template <typename... Cols>
        class streamed_inmemory_table : public input_table_implementation<Cols...>
        {
        public:
            using row = typename table_implementation<Cols...>::row;

            streamed_inmemory_table(const inmemory_table<Cols...>& p_tbl)
                : m_begin(p_tbl.begin()), m_end(p_tbl.end())
            {
            }

            bool next(row& p_row)
            {
                if (m_begin == m_end)
                {
                    return false;
                }
                p_row = *m_begin;
                ++m_begin;
                return true;
            }

        private:
            typename inmemory_table<Cols...>::const_iterator m_begin;
            typename inmemory_table<Cols...>::const_iterator m_end;
        };

        template <typename... Cols>
        class streamed_outmemory_table : public output_table_implementation<Cols...>
        {
        public:
            using row = typename table_implementation<Cols...>::row;

            streamed_outmemory_table(inmemory_table<Cols...>& p_tbl)
                : m_tbl(p_tbl)
            {
            }

            streamed_outmemory_table& operator<<(const row& p_row)
            {
                m_tbl.push_back(p_row);
                return *this;
            }

        private:
            inmemory_table<Cols...>& m_tbl;
            typename inmemory_table<Cols...>::const_iterator m_end;
        };

    }
    // namespace detail

    struct table_utils
    {
        template <std::size_t I0, std::size_t I1, std::size_t O = 0,
                  typename... Ts, template <typename...> class T,
                  typename... Us, template <typename...> class U>
        static 
        typename std::enable_if<I0 == I1, void>::type
        copy(const T<Ts...>& p_src, U<Us...>& p_dest)
        {
        }

        template <std::size_t I0, std::size_t I1, int O = 0,
                  typename... Ts, template <typename...> class T,
                  typename... Us, template <typename...> class U>
        static 
        typename std::enable_if<I0 < I1, void>::type
        copy(const T<Ts...>& p_src, U<Us...>& p_dest)
        {
            static_assert(I0+O >= 0);
            static_assert(std::is_base_of<std::tuple<Ts...>, T<Ts...>>::value);
            static_assert(std::is_base_of<std::tuple<Us...>, U<Us...>>::value);
            static_assert(I0 < std::tuple_size<T<Ts...>>::value);
            static_assert(I1 <= std::tuple_size<T<Ts...>>::value);
            static_assert(I0+O < std::tuple_size<U<Us...>>::value);
            static_assert(I1+O <= std::tuple_size<U<Us...>>::value);

            std::get<I0+O>(p_dest) = std::get<I0>(p_src);
            copy<I0+1, I1, O>(p_src, p_dest);
        }

        template <typename... Ts, template <typename...> class T, template <typename...> class U, template <typename...> class V>
        static void merge(T<Ts...>& p_lhs, U<Ts...>& p_rhs, V<Ts...>& p_res,
                          std::function<bool(const std::tuple<Ts...>&, const std::tuple<Ts...>&)> p_less,
                          std::function<void(const std::tuple<Ts...>&, const std::tuple<Ts...>&, std::tuple<Ts...>&)> p_combine)
        {
            static_assert(std::is_base_of<detail::input_table_implementation<Ts...>, T<Ts...>>::value);
            static_assert(std::is_base_of<detail::input_table_implementation<Ts...>, U<Ts...>>::value);
            static_assert(std::is_base_of<detail::output_table_implementation<Ts...>, V<Ts...>>::value);

            using row_type = std::tuple<Ts...>;

            row_type x;
            bool x_more = p_lhs.next(x);
            row_type y;
            bool y_more = p_rhs.next(y);

            row_type z;
            while (x_more && y_more)
            {
                if (p_less(x, y))
                {
                    p_res << x;
                    x_more = p_lhs.next(x);
                    continue;
                }
                if (p_less(y, x))
                {
                    p_res << y;
                    y_more = p_rhs.next(y);
                    continue;
                }
                p_combine(x, y, z);
                p_res << z;
                x_more = p_lhs.next(x);
                y_more = p_rhs.next(y);
            }
            while (x_more)
            {
                p_res << x;
                x_more = p_lhs.next(x);
            }
            while (y_more)
            {
                p_res << y;
                y_more = p_rhs.next(y);
            }
        }

        template <typename... Ts, template <typename...> class T>
        static void for_each(T<Ts...>& p_src, std::function<void(const std::tuple<Ts...>&)> p_func)
        {
            static_assert(std::is_base_of<detail::input_table_implementation<Ts...>, T<Ts...>>::value);

            using row_in = std::tuple<Ts...>;

            row_in x;

            while(p_src.next(x))
            {
                p_func(x);
            }
        }

        template <typename... Ts, template <typename...> class T,
                  typename... Us, template <typename...> class U>
        static void for_each_2(T<Ts...>& p_lhs, U<Us...>& p_rhs, std::function<void(const std::tuple<Ts...>&,
                                                                                    const std::tuple<Us...>&)> p_func)
        {
            static_assert(std::is_base_of<detail::input_table_implementation<Ts...>, T<Ts...>>::value);
            static_assert(std::is_base_of<detail::input_table_implementation<Us...>, U<Us...>>::value);

            using lhs_row = std::tuple<Ts...>;
            using rhs_row = std::tuple<Us...>;

            lhs_row x;
            rhs_row y;

            bool more_x = p_lhs.next(x);
            bool more_y = p_rhs.next(y);

            while(more_x && more_y)
            {
                p_func(x, y);
                more_x = p_lhs.next(x);
                more_y = p_rhs.next(y);
            }
            if (more_x || more_y)
            {
                throw std::runtime_error("for_each: iterators have different lengths");
            }
        }

        template <typename... Ts, template <typename...> class T,
                  typename... Us, template <typename...> class U>
        static void map(std::function<void(const std::tuple<Ts...>&,std::tuple<Us...>&)> p_func, T<Ts...>& p_src, U<Us...>& p_dest)
        {
            static_assert(std::is_base_of<detail::input_table_implementation<Ts...>, T<Ts...>>::value);
            static_assert(std::is_base_of<detail::output_table_implementation<Us...>, U<Us...>>::value);

            using row_in = std::tuple<Ts...>;
            using row_out = std::tuple<Us...>;

            row_in x;
            row_out y;

            while(p_src.next(x))
            {
                p_func(x, y);
                p_dest << y;
            }
        }

    };

    template <typename... Ts>
    class table
    {
    public:
        using tuple_type = std::tuple<Ts...>;

        using basic_istream_table = detail::streamed_input_table<Ts...>;

        using basic_ostream_table = detail::streamed_output_table<Ts...>;

        using basic_inmemory_table = detail::inmemory_table<Ts...>;

        using basic_read_iterator = detail::streamed_inmemory_table<Ts...>;

        using basic_write_iterator = detail::streamed_outmemory_table<Ts...>;

        static void write(std::ostream& p_out, const basic_inmemory_table& p_table)
        {
            detail::streamed_output_table<Ts...> out(p_out);
            for (auto itr = p_table.begin(); itr != p_table.end(); ++itr)
            {
                out << *itr;
            }
        }

        static void write(std::ostream& p_out, const basic_inmemory_table& p_table, std::initializer_list<std::string> p_labels)
        {
            detail::streamed_output_table<Ts...> out(p_out, p_labels);
            for (auto itr = p_table.begin(); itr != p_table.end(); ++itr)
            {
                out << *itr;
            }
        }

    };

}
// namespace varoom

#endif // VAROOM_UTIL_TABLE_HPP
