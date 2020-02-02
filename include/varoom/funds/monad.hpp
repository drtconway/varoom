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
            static constexpr bool is_instance = false;

            template <typename U, typename T, typename X>
            static F<U> fmap(X p_m, F<T> p_x);
        };

        template <template<typename> class F, typename T, typename X>
        auto fmap(X x, F<T> t) -> F<decltype(x(*reinterpret_cast<const T*>(NULL)))>
        {
            static_assert(functor<F>::is_instance);
            using U = decltype(x(*reinterpret_cast<const T*>(NULL)));
            static_assert(std::is_convertible<X, std::function<U(T)>>::value,
                          "applicative requires a function type U(T)");
            F<U> res = functor<F>::template fmap<U>(x, t);
            return res;
        }

        template <template<typename> class A>
        struct applicative
        {
            static constexpr bool is_instance = false;

            //template <typename T>
            //static A<T> pure(T p_x);

            //template <typename U, typename T, typename X>
            //static A<U> apply(A<X> p_m, A<T> p_x);
        };

        template <template<typename> class A, typename T>
        A<T> pure(const T& t)
        {
            static_assert(applicative<A>::is_instance);
            return applicative<A>::template pure<T>(t);
        }

        template <template<typename> class A, typename T, typename X>
        auto apply(A<X> x, A<T> t) -> A<decltype((*reinterpret_cast<const X*>(NULL))(*reinterpret_cast<const T*>(NULL)))>
        {
            static_assert(applicative<A>::is_instance);
            using U = decltype((*reinterpret_cast<const X*>(NULL))(*reinterpret_cast<const T*>(NULL)));
            static_assert(std::is_convertible<X, std::function<U(T)>>::value,
                          "applicative requires a function type U(T)");
            A<U> res = applicative<A>::template apply<U>(x, t);
            return res;
        }

        template <template<typename> class M>
        struct monoid
        {
            static constexpr bool is_instance = false;

            //template <typename T>
            //static M<T> empty();

            //template <typename T>
            //static M<T> append(M<T>, M<T>);
        };

        template <template<typename> class M, typename T>
        M<T> empty()
        {
            static_assert(monoid<M>::is_instance);
            return monoid<M>::template empty<T>();
        }

        template <template<typename> class M, typename T>
        M<T> append(M<T> x, M<T> y)
        {
            static_assert(monoid<M>::is_instance);
            return monoid<M>::template append<T>(x, y);
        }

        template <template<typename> class M>
        struct monad
        {
            static constexpr bool is_instance = false;

            //template <typename T>
            //static M<T> yield(T p_x)

            //template <typename U, typename T, typename X>
            //static M<U> for_each(X p_m, M<T> p_x)
        };

        template <template<typename> class M, typename T>
        M<T> yield(T t)
        {
            static_assert(monad<M>::is_instance);
            return monad<M>::template yield<T>(t);
        }

        namespace detail
        {
            template <template<typename> class M, typename T>
            const T* unwrap_type(const M<T> mt)
            {
                return NULL;
            }

            template <typename T, typename U>
            struct same_functor : std::false_type {};

            template <template<typename> class F, typename T, typename U>
            struct same_functor<F<T>,F<U>> : std::true_type
            {
                using lhs_type = T;
                using rhs_type = U;
            };
        }
        // namespace detail

        template <template<typename> class M, typename T, typename X>
        auto bind(X x, M<T> t) -> decltype((*reinterpret_cast<const X*>(NULL))(*reinterpret_cast<const T*>(NULL)))
        {
            static_assert(monad<M>::is_instance);
            using MU = decltype(x(*reinterpret_cast<T*>(NULL)));
            using V = detail::same_functor<M<T>,MU>;
            using U = typename V::rhs_type;
            static_assert(std::is_convertible<X, std::function<M<U>(T)>>::value,
                          "applicative requires a function type U(T)");
            M<U> res = monad<M>::template bind<U>(x, t);
            return res;
        }

        template <template<typename> class M, typename T, typename X>
        typename std::enable_if<std::is_convertible<X, std::function<void(T)>>::value, void>::type
        for_each(X x, M<T> t)
        {
            static_assert(monad<M>::is_instance);
            monad<M>::template for_each(x, t);
        }
    }
    // namespace funds
}
// namespace varoom

#endif // VAROOM_FUNDS_MONAD_HPP
