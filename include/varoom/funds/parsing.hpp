#ifndef VAROOM_FUNDS_PARSING_HPP
#define VAROOM_FUNDS_PARSING_HPP

#ifndef VAROOM_FUNDS_LIST_HPP
#include "varoom/funds/list.hpp"
#endif

#ifndef VAROOM_FUNDS_maybe_HPP
#include "varoom/funds/maybe.hpp"
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
        struct functor<parser>
        {
            static constexpr bool is_instance = true;

            template <typename U, typename T, typename X>
            static parser<U> fmap(X x, parser<T> p)
            {
                static_assert(std::is_convertible<X, std::function<U(T)>>::value);
                return [=](symbols s) {
                    pair_list<T,symbols> ts = parse(p, s);
                    return functor<list>::template fmap<std::pair<U,symbols>>([=](std::pair<T,symbols> t) {
                        std::pair<U,symbols> us(x(t.first), t.second);
                        return us;
                    }, ts);
                };
            }
        };
        
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

        template <typename T>
        parser<T> appendx(parser<T> x, parser<T> y)
        {
            return [=](symbols s) {
                pair_list<T,symbols> x_res = parse(x, s);
                if (!x_res.empty())
                {
                    return x_res;
                }
                else
                {
                    return parse(y, s);
                }
            };
        }

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

        template <typename T, typename X>
        auto operator>>=(parser<T> p, X f) -> decltype(f(*reinterpret_cast<const T*>(NULL)))
        {
            using MU = decltype(f(*reinterpret_cast<T*>(NULL)));
            using V = detail::same_functor<parser<T>,MU>;
            using U = typename V::rhs_type;
            return monad<parser>::template bind<U>(f, p);
        }

        template <typename T, typename X>
        parser<T> sat(parser<T> p, X x)
        {
            static_assert(std::is_convertible<X, std::function<bool(T)>>::value);
            return (p >>= [=](T t) {
                if (x(t))
                {
                    return yield<parser,T>(t);
                }
                else
                {
                    return empty<parser,T>();
                }
            });
        }

        template <typename T>
        parser<list<T>> many0(parser<T> p)
        {
            return appendx((p >>= [=](T t) {
                return (many0(p) >>= [=](list<T> ts) {
                    return yield<parser,list<T>>(list<T>(t, ts));
                });
            }), yield<parser,list<T>>(list<T>{}));
        }

        template <typename T>
        parser<list<T>> many1(parser<T> p)
        {
            return (p >>= [=](T t) {
                return (many0(p) >>= [=](list<T> ts) {
                    return yield<parser,list<T>>(list<T>(t, ts));
                });
            });
        }

        template <typename T>
        parser<maybe<T>> optional(parser<T> p)
        {
            return appendx(p >>= [](T t) {
                return varoom::funds::yield<parser,maybe<T>>(maybe<T>(t));
            }, varoom::funds::yield<parser,maybe<T>>(maybe<T>{}));
        }

        template <typename T>
        parser<T> alt(std::initializer_list<parser<T>> p_alts)
        {
            if (p_alts.begin() == p_alts.end())
            {
                return empty<parser,T>();
            }
            auto itr = p_alts.begin();
            parser<T> p = *itr;
            while (++itr != p_alts.end())
            {
                p = append(p, *itr);
            }
            return p;
        }

        template <typename T>
        parser<T> altx(std::initializer_list<parser<T>> p_alts)
        {
            if (p_alts.begin() == p_alts.end())
            {
                return empty<parser,T>();
            }
            auto itr = p_alts.begin();
            parser<T> p = *itr;
            while (++itr != p_alts.end())
            {
                p = appendx(p, *itr);
            }
            return p;
        }

        parser<char> sym()
        {
            return [](symbols s) {
                auto itr = s.begin();
                if (itr == s.end())
                {
                    return pair_list<char,symbols>();
                }
                else
                {
                    char t = *itr;
                    symbols r = symbols(++itr, s.end());
                    return pair_list<char,symbols>(std::make_pair(t, r));
                }
            };
        }

        parser<char> sym(char c)
        {
            return sat(sym(), [=](char d) { return c == d; });
        }

        parser<char> oneof(std::string p_chars)
        {
            return sat(sym(), [=](char c) {
                for (size_t i = 0; i < p_chars.size(); ++i)
                {
                    if (c == p_chars[i])
                    {
                        return true;
                    }
                }
                return false;
            });
        }

        std::string list_to_string(list<char> ts)
        {
            std::string s;
            while (!ts.empty())
            {
                s.push_back(ts.head());
                ts = ts.tail();
            }
            return s;
        }
    }
    // namespace funds
}
// namespace varoom

#endif // VAROOM_FUNDS_PARSING_HPP
