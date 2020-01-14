#ifndef VAROOM_FUNDS_LIST_HPP
#define VAROOM_FUNDS_LIST_HPP

#include <memory>

#ifndef VAROOM_FUNDS_MONAD_HPP
#include "varoom/funds/monad.hpp"
#endif

namespace varoom
{
    namespace funds
    {
        template <typename T>
        class list
        {
        private:
            struct item
            {
                T value;
                std::shared_ptr<item> next;

                item(T p_value)
                    : value(p_value)
                {
                }

                item(T p_value, const std::shared_ptr<item>& p_next)
                    : value(p_value), next(p_next)
                {
                }
            };
            using item_ptr = std::shared_ptr<item>;

            list(const item_ptr& p_list)
                : m_head(p_list)
            {
            }

        public:
            list() {}
            
            list(T p_x)
                : m_head(std::make_shared<item>(p_x, item_ptr()))
            {
            }

            list(T p_x, const list& p_tail)
                : m_head(std::make_shared<item>(p_x, p_tail.m_head))
            {
            }

            bool empty() const
            {
                return m_head == nullptr;
            }

            T head() const
            {
                return m_head->value;
            }

            list tail() const
            {
                return list(m_head->next);
            }

            static list concat(const list& p_lhs, const list& p_rhs)
            {
                if (p_lhs.empty())
                {
                    return p_rhs;
                }
                if (p_rhs.empty())
                {
                    return p_lhs;
                }
                return list(p_lhs.head(), concat(p_lhs.tail(), p_rhs));
            }

            static list flatten(list<list<T>> xss)
            {
                list<T> xs;
                while (!xss.empty())
                {
                    xs = concat(xs, xss.head());
                    xss = xss.tail();
                }
                return xs;
            }

        private:
            item_ptr m_head;
        };

        template <>
        struct functor<list>
        {
            using is_implemented = std::true_type;

            template <typename U, typename T, typename X>
            static list<U> fmap(X p_m, list<T> p_xs)
            {
                static_assert(std::is_convertible<X, std::function<U(T)>>::value, 
                              "fmap requires a function type U(T)");

                if (p_xs.empty())
                {
                    return list<U>{};
                }
                else
                {
                    list<U> tl = fmap<U,T,X>(p_m, p_xs.tail());
                    return list<U>(p_m(p_xs.head()), tl);
                }
            }
        };

        template <typename U, typename T, typename X>
        list<U> fmap(X p_xform, list<T> p_xs)
        {
            return functor<list>::fmap<U,T,X>(p_xform, p_xs);
        }

        template <typename T, typename X>
        list<T> filter(X p_pred, list<T> p_xs)
        {
            static_assert(std::is_convertible<X, std::function<bool(T)>>::value, 
                          "filter requires a function type bool(T)");

            if (p_xs.empty())
            {
                return p_xs;
            }
            if (p_pred(p_xs.head()))
            {
                return list<T>(p_xs.head(), filter(p_pred, p_xs.tail()));
            }
            else
            {
                return filter(p_pred, p_xs.tail());
            }
        }

        template <>
        struct applicative<list>
        {
            using is_implemented = std::true_type;

            template <typename T>
            static list<T> pure(T p_x)
            {
                return list<T>(p_x);
            }

            template <typename T, typename U, typename X>
            static list<U> apply(list<X> p_ms, list<T> p_xs)
            {
                static_assert(std::is_convertible<X, std::function<U(T)>>::value, 
                              "apply requires a function type U(T)");

                if (p_ms.empty())
                {
                    return list<U>{};
                }
                else
                {
                    return list<U>::concat(functor<list>::fmap<U,T,X>(p_ms.head(), p_xs), apply<T,U,X>(p_ms.tail(), p_xs));
                }
            }
        };

        template <typename T>
        list<T> pure(T p_x)
        {
            return applicative<list>::pure<T>(p_x);
        }

        template <typename U, typename T, typename X>
        list<U> apply(list<X> p_xs, list<T> p_ts)
        {
            return applicative<list>::apply<T,U,X>(p_xs, p_ts);
        }

        template <>
        struct monad<list>
        {
            using is_implemented = std::true_type;

            template <typename T>
            static list<T> yield(T p_x)
            {
                return list<T>(p_x);
            }

            template <typename U, typename T, typename X>
            static list<U> for_each(X p_x, list<T> p_ts)
            {
                static_assert(std::is_convertible<X, std::function<list<U>(T)>>::value, 
                              "for_each requires a function type U(T)");

                list<list<U>> uss = functor<list>::fmap<U>(p_x, p_ts);
                return list<U>::flatten(uss);
            }

            template <typename T, typename X>
            static void for_each(X x, list<T> ts)
            {
                static_assert(std::is_convertible<X, std::function<void(T)>>::value, 
                              "for_each requires a function type void(T)");

                while (!ts.empty())
                {
                    x(ts.head());
                    ts = ts.tail();
                }
            }
        };

        template <typename T>
        list<T> yield(T x)
        {
            return monad<list>::yield(x);
        }

        template <typename U, typename T, typename X>
        static list<U> for_each(X p_x, list<T> p_ts)
        {
            return monad<list>::for_each<U,T,X>(p_x, p_ts);
        }

        template <typename T, typename X>
        static void for_each(X p_x, list<T> p_ts)
        {
            monad<list>::for_each<T,X>(p_x, p_ts);
        }
    }
    // namespace funds
}
// namespace varoom

#endif // VAROOM_FUNDS_LIST_HPP
