#ifndef VAROOM_FUNDS_MONAD_HPP
#define VAROOM_FUNDS_MONAD_HPP

#include <functional>

namespace varoom
{
    namespace funds
    {
        template <template<typename> class F>
        struct functor
        {
            using is_implemented = std::false_type;

            template <typename T, typename U, typename X>
            static F<U> fmap(X p_m, F<T> p_x)
            {
                static_assert(std::is_convertible<X, std::function<U(T)>>::value, 
                              "fmap requires a function type U(T)");
            }
        };

        template <template<typename> class A>
        struct applicative
        {
            using is_implemented = std::false_type;

            template <typename T>
            static A<T> pure(T p_x)
            {
            }

            template <typename T, typename U, typename X>
            static A<U> apply(A<X> p_m, A<T> p_x)
            {
                static_assert(std::is_convertible<X, std::function<U(T)>>::value, 
                              "applicative requires a function type U(T)");
            }
        };

        template <template<typename> class M>
        struct monad
        {
            using is_implemented = std::false_type;

            template <typename T>
            static M<T> yield(T p_x)
            {
            }

            template <typename T, typename U, typename X>
            static M<U> for_each(X p_m, M<T> p_x)
            {
                static_assert(std::is_convertible<X, std::function<M<U>(T)>>::value, 
                              "for_each requires a function type M<U>(T)");
            }
        };
    }
    // namespace funds
}
// namespace varoom

#endif // VAROOM_FUNDS_MONAD_HPP
