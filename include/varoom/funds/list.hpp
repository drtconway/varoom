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

                item(T p_value, std::shared_ptr<item>& p_next)
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
                : m_head(std::make_shared(p_x));
            {
            }

            list(T p_x, const list& p_tail)
                : m_head(std::make_shared(p_x, p_tail.m_head))
            {
            }

            bool empty() const
            {
                return m_head;
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

        private:
            item_ptr m_head;
        };

        template <>
        struct functor<list>
        {
            using is_implemented = std::true_type;

            template <typename T, typename U, typename M>
            static list<U> fmap(M p_m, list<T> p_xs)
            {
                static_assert(std::is_convertible<M, std::function<U(T)>>::value, 
                              "fmap requires a function type U(T)");

                if (p_xs.empty())
                {
                    return list<U>{};
                }
                else
                {
                    return list<U>{p_m(p_xs.head), fmap(p_m, p_xs.tail())};
                }
            }
        };

        template <>
        struct applicative<list>
        {
            using is_implemented = std::true_type;

            template <typename T>
            static list<T> pure(T p_x)
            {
                return list(p_x);
            }

            template <typename T, typename U, typename M>
            static list<U> apply(list<M> p_ms, list<T> p_xs)
            {
                static_assert(std::is_convertible<M, std::function<U(T)>>::value, 
                              "fmap requires a function type U(T)");

                if (p_ms.empty())
                {
                    return list<U>{};
                }
                else
                {
                    retrun list<U>::concat(functor<A>::fmap(p_ms.head(), p_xs), apply(p_ms.tail(), p_xs));
                }
            }
        };
    }
    // namespace funds
}
// namespace varoom

#endif VAROOM_FUNDS_LIST_HPP
