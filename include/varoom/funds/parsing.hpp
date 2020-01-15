#ifndef VAROOM_FUNDS_PARSING_HPP
#define VAROOM_FUNDS_PARSING_HPP

#ifndef VAROOM_FUNDS_LIST_HPP
#include "varoom/funds/list.hpp"
#endif

namespace varoom
{
    namespace funds
    {
        using state = std::string;

        template <typename T, typename U>
        using pair_list = varoom::funds::list<std::pair<T,U>>;

        class symbols
        {
        public:
            using const_iterator = std::string::const_iterator;

            symbols(const std::string& p_str)
                : m_begin(p_str.begin()), m_end(p_str.end())
            {
            }

            symbols(const_iterator p_begin, const_iterator p_end)
                : m_begin(p_begin), m_end(p_end)
            {
            }

            const_iterator begin() const
            {
                return m_begin;
            }

            const_iterator end() const
            {
                return m_end;
            }

        private:
            const_iterator m_begin;
            const_iterator m_end;
        };

        template <typename T>
        struct parser : std::function<pair_list<T,symbols>(symbols)>
        {
            template <typename X>
            parser(X x) : std::function<pair_list<T,symbols>(symbols)>(x) {}
        };

        template <typename T>
        pair_list<T,symbols> parse(parser<T> p, symbols s)
        {
            return p(s);
        }
        
        template <>
        struct monoid<parser>
        {
            static constexpr bool is_instance = true;

            template <typename T>
            static parser<T> empty()
            {
                return [](symbols s) {
                    return pair_list<T,symbols>();
                };
            }

            template <typename T>
            static parser<T> append(parser<T> lhs, parser<T> rhs)
            {
                return [=](symbols s) {
                    pair_list<T,symbols> lhs_res = parse(lhs, s);
                    pair_list<T,symbols> rhs_res = parse(rhs, s);
                    return pair_list<T,symbols>::concat(lhs_res, rhs_res);
                };
            }
        };

        template<>
        struct monad<parser>
        {
            static constexpr bool is_instance = true;

            template <typename T>
            static parser<T> yield(T t)
            {
                return [=](symbols s) {
                    return pair_list<T,symbols>(std::make_pair(t, s));
                };
            }

            template<typename U, typename T, typename X>
            static parser<U> bind(X x, parser<T> t)
            {
                static_assert(std::is_convertible<X, std::function<parser<U>(T)>>::value);

                return [=](symbols s) {
                    list<pair_list<U,symbols>> uss = fmap([=](std::pair<T,symbols> p) {
                        return parse(x(p.first), p.second);
                    }, parse(t, s));
                    return pair_list<U,symbols>::flatten(uss);
                };
            }
        };

    }
    // namespace funds
}
// namespace varoom

#endif // VAROOM_FUNDS_PARSING_HPP
