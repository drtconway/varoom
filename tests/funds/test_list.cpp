#include "varoom/funds/list.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE funds_list tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
    template <typename T>
    size_t length(varoom::funds::list<T> xs)
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
    T access(varoom::funds::list<T> xs, size_t i)
    {
        while (i > 0)
        {
            xs = xs.tail();
            i -= 1;
        }
        return xs.head();
    }

    template<typename T>
    varoom::funds::list<T> make(std::initializer_list<T> xs)
    {
        std::vector<T> v(xs.begin(), xs.end());

        varoom::funds::list<T> ys;
        for (auto itr = v.rbegin(); itr != v.rend(); ++itr)
        {
            ys = varoom::funds::list<T>(*itr, ys);
        }
        return ys;
    }
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( cons1 )
{
    using L = varoom::funds::list<int>;

    L e;

    BOOST_CHECK_EQUAL(e.empty(), true);
}

BOOST_AUTO_TEST_CASE( cons2 )
{
    using L = varoom::funds::list<int>;

    L s{1};

    BOOST_CHECK_EQUAL(s.empty(), false);
    BOOST_CHECK_EQUAL(s.head(), 1);
    BOOST_CHECK_EQUAL(s.tail().empty(), true);
    BOOST_CHECK_EQUAL(length(s), 1);
}

BOOST_AUTO_TEST_CASE( fmap1 )
{
    using L = varoom::funds::list<int>;
    using M = varoom::funds::list<double>;

    L s{1, L(2, L(3))};

    BOOST_CHECK_EQUAL(length(s), 3);
    BOOST_CHECK_EQUAL(access(s, 0), 1);
    BOOST_CHECK_EQUAL(access(s, 1), 2);
    BOOST_CHECK_EQUAL(access(s, 2), 3);

    M t = varoom::funds::fmap([=](int p_x) {
        return 1.0/double(p_x);
    }, s);

    BOOST_CHECK_EQUAL(length(t), 3);
    BOOST_CHECK_EQUAL(access(t, 0), 1.0/1.0);
    BOOST_CHECK_EQUAL(access(t, 1), 1.0/2.0);
    BOOST_CHECK_EQUAL(access(t, 2), 1.0/3.0);
}

BOOST_AUTO_TEST_CASE( pure1 )
{
    using L = varoom::funds::list<int>;

    L s = varoom::funds::pure<varoom::funds::list>(1);

    BOOST_CHECK_EQUAL(length(s), 1);
    BOOST_CHECK_EQUAL(access(s, 0), 1);
}

BOOST_AUTO_TEST_CASE( apply1 )
{
    using L = varoom::funds::list<int>;
    using M = varoom::funds::list<std::function<double(int)>>;
    using N = varoom::funds::list<double>;

    L s = varoom::funds::pure<varoom::funds::list>(2);
    M t = varoom::funds::pure<varoom::funds::list>(std::function<double(int)>{[=](int x) { return 1.0/double(x); }});
    N u = varoom::funds::apply(t, s);

    BOOST_CHECK_EQUAL(length(u), 1);
    BOOST_CHECK_EQUAL(access(u, 0), 0.5);
}

BOOST_AUTO_TEST_CASE( apply2 )
{
    using L = varoom::funds::list<int>;
    using M = varoom::funds::list<std::function<double(int)>>;
    using N = varoom::funds::list<double>;

    L s = make({1, 2, 3});
    M t = make({std::function<double(int)>([=](int x) { return 1.0/double(x); }),
                std::function<double(int)>([=](int x) { return double(x) * double(x); })});
    N u = varoom::funds::apply(t, s);

    BOOST_CHECK_EQUAL(length(u), 6);
    BOOST_CHECK_EQUAL(access(u, 0), 1.0/1.0);
    BOOST_CHECK_EQUAL(access(u, 1), 1.0/2.0);
    BOOST_CHECK_EQUAL(access(u, 2), 1.0/3.0);
    BOOST_CHECK_EQUAL(access(u, 3), 1.0*1.0);
    BOOST_CHECK_EQUAL(access(u, 4), 2.0*2.0);
    BOOST_CHECK_EQUAL(access(u, 5), 3.0*3.0);
}

BOOST_AUTO_TEST_CASE( monad1 )
{
    using L = varoom::funds::list<int>;

    L s = varoom::funds::yield<varoom::funds::list>(2);

    BOOST_CHECK_EQUAL(length(s), 1);
    BOOST_CHECK_EQUAL(access(s, 0), 2);

    std::vector<int> xs;
    varoom::funds::for_each([&](int x) mutable {
        xs.push_back(x);
    }, s);

    BOOST_CHECK_EQUAL(xs.size(), 1);
    BOOST_CHECK_EQUAL(xs[0], 2);
}

BOOST_AUTO_TEST_CASE( monad2 )
{
    using L = varoom::funds::list<int>;
    using M = varoom::funds::list<double>;

    static_assert(varoom::funds::detail::same_functor<L,M>::value);
    static_assert(std::is_same<varoom::funds::detail::same_functor<L,M>::lhs_type,int>::value);
    static_assert(std::is_same<varoom::funds::detail::same_functor<L,M>::rhs_type,double>::value);

    L s = varoom::funds::yield<varoom::funds::list>(2);

    BOOST_CHECK_EQUAL(length(s), 1);
    BOOST_CHECK_EQUAL(access(s, 0), 2);

    std::function<M(int)> f = [&](int x) {
        return varoom::funds::yield<varoom::funds::list>(double(x));
    };

    M t = varoom::funds::bind(f, s);
    std::vector<double> xs;
    varoom::funds::for_each([&](double x) mutable {
        xs.push_back(x);
    }, t);

    BOOST_CHECK_EQUAL(xs.size(), 1);
    BOOST_CHECK_EQUAL(xs[0], 2.0);
}

