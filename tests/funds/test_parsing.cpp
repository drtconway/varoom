#include "varoom/funds/parsing.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE funds_parsing tests
#include <boost/test/unit_test.hpp>

using namespace varoom::funds;

namespace // anonymous
{
    template <typename T>
    size_t length(list<T> xs)
    {
        size_t n = 0;
        while (!xs.empty())
        {
            xs = xs.tail();
            n += 1;
        }
        return n;
    }

    template <typename T>
    T access(list<T> xs, size_t i)
    {
        while (i > 0)
        {
            xs = xs.tail();
            i -= 1;
        }
        return xs.head();
    }

    template <typename T, typename X>
    auto operator>>=(parser<T> p, X f) -> decltype(f(*reinterpret_cast<const T*>(NULL)))
    {
        using MU = decltype(f(*reinterpret_cast<T*>(NULL)));
        using V = detail::same_functor<parser<T>,MU>;
        using U = typename V::rhs_type;
        return monad<parser>::template bind<U>(f, p);
    }
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( cons1 )
{
    using P = parser<int>;
    P p = empty<parser,int>();
    P q = append(p, p);

    auto r = parse(q, symbols(""));
    BOOST_CHECK_EQUAL(length(r), 0);
}

BOOST_AUTO_TEST_CASE( cons2 )
{
    using P = parser<int>;
    using Q = parser<double>;
    P p = yield<parser,int>(2);
    Q q = (p >>= [](int x) {
            return yield<parser,double>(3.14*x*x);
          });

    auto r = parse(q, symbols(""));
    BOOST_CHECK_EQUAL(length(r), 1);
    BOOST_CHECK_EQUAL(access(r, 0).first, 3.14*2*2);
}
