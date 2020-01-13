#ifndef VAROOM_FUNDS_MONAD_HPP
#define VAROOM_FUNDS_MONAD_HPP

namespace varoom
{
    namespace funds
    {
        template <template<typename> F>
        struct functor
        {
            using is_implemented = std::false_type;

            template <typename T, typename U, typename X>
            static F<U> fmap(X p_m, F<T> p_x)
            {
                static_assert(std::is_convertible<X, std::function<U(T)>>::value, 
                              "fmap requires a function type U(T)");
                static_assert(std::false_type, "fmap is not implemented");
            }
        };

        template <template<typename> A>
        struct applicative
        {
            using is_implemented = std::false_type;

            template <typename T>
            static A<T> pure(T p_x)
            {
                static_assert(std::false_type, "pure is not implemented");
            }

            template <typename T, typename U, typename X>
            static A<U> apply(A<X> p_m, A<T> p_x)
            {
                static_assert(std::is_convertible<X, std::function<U(T)>>::value, 
                              "fmap requires a function type U(T)");
                static_assert(std::false_type, "pure is not implemented");
            }
        };
    }
    // namespace funds
}
// namespace varoom

#endif // VAROOM_FUNDS_MONAD_HPP
