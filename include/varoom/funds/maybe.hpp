#ifndef VAROOM_FUNDS_MAYBE_HPP
#define VAROOM_FUNDS_MAYBE_HPP

#include <memory>

#ifndef VAROOM_FUNDS_MONAD_HPP
#include "varoom/funds/monad.hpp"
#endif

namespace varoom
{
    namespace funds
    {
        template <typename T>
        class maybe
        {
        private:
            using item = T;
            using item_ptr = std::shared_ptr<item>;

        public:
            maybe() {}
            
            maybe(T p_x)
                : m_item(std::make_shared<item>(p_x))
            {
            }

            bool nothing() const
            {
                return m_item == nullptr;
            }

            T just() const
            {
                return m_item->value;
            }

        private:
            item_ptr m_item;
        };

        template <>
        struct functor<maybe>
        {
            static constexpr bool is_instance = true;

            template <typename U, typename T, typename X>
            static maybe<U> fmap(X p_m, maybe<T> p_xs)
            {
                static_assert(std::is_convertible<X, std::function<U(T)>>::value, 
                              "fmap requires a function type U(T)");

                if (p_xs.nothing())
                {
                    return maybe<U>{};
                }
                else
                {
                    return maybe<U>(p_m(p_xs.just()));
                }
            }
        };

        template <>
        struct monad<maybe>
        {
            static constexpr bool is_instance = true;

            template <typename T>
            static maybe<T> yield(T p_x)
            {
                return maybe<T>(p_x);
            }

            template <typename U, typename T, typename X>
            static maybe<U> bind(X p_x, maybe<T> p_ts)
            {
                static_assert(std::is_convertible<X, std::function<maybe<U>(T)>>::value, 
                              "for_each requires a function type U(T)");

                if (p_ts.nothing())
                {
                    return monad<maybe>::template yield<U>(maybe<U>{});
                }
                else
                {
                    return monad<maybe>::template yield<U>(maybe<U>{x(p_ts.just())});
                }
            }

            template <typename T, typename X>
            static void for_each(X x, list<T> ts)
            {
                static_assert(std::is_convertible<X, std::function<void(T)>>::value, 
                              "for_each requires a function type void(T)");

                if (!ts.empty())
                {
                    x(ts.just());
                }
            }
        };
    }
    // namespace funds
}
// namespace varoom

#endif // VAROOM_FUNDS_MAYBE_HPP
