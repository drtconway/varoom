#include "varoom/funds/list.hpp"
#include "varoom/funds/parsing.hpp"

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

    bool oneof(const std::string& p_str, char c)
    {
        for (auto itr = p_str.begin(); itr != p_str.end(); ++itr)
        {
            if (*itr == c)
            {
                return true;
            }
        }
        return false;
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

    M t = varoom::funds::fmap<double>([=](int p_x) {
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

    L s = varoom::funds::pure(1);

    BOOST_CHECK_EQUAL(length(s), 1);
    BOOST_CHECK_EQUAL(access(s, 0), 1);
}

BOOST_AUTO_TEST_CASE( apply1 )
{
    using L = varoom::funds::list<int>;
    using M = varoom::funds::list<std::function<double(int)>>;
    using N = varoom::funds::list<double>;

    L s = varoom::funds::pure(2);
    M t = varoom::funds::pure(std::function<double(int)>{[=](int x) { return 1.0/double(x); }});
    N u = varoom::funds::apply<double>(t, s);

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
    N u = varoom::funds::apply<double>(t, s);

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

    L s = varoom::funds::yield(2);

    BOOST_CHECK_EQUAL(length(s), 1);
    BOOST_CHECK_EQUAL(access(s, 0), 2);

    std::vector<int> xs;
    varoom::funds::for_each([&](int x) mutable {
        xs.push_back(x);
    }, s);

    BOOST_CHECK_EQUAL(xs.size(), 1);
    BOOST_CHECK_EQUAL(xs[0], 2);
}

BOOST_AUTO_TEST_CASE( parsing1 )
{
    using T = varoom::funds::parser<std::string,int>;

    std::string s = "wibble";
    T t = varoom::funds::yield<int,std::string>(1);
    auto u = varoom::funds::parse(t, s);
    BOOST_CHECK_EQUAL(length(u), 1);
    BOOST_CHECK_EQUAL(access(u, 0).first, 1);
    BOOST_CHECK_EQUAL(access(u, 0).second, s);
}

BOOST_AUTO_TEST_CASE( parsing2 )
{
    using T = varoom::funds::parser<std::string,int>;
    using U = varoom::funds::parser<std::string,double>;

    std::string s = "wibble";
    T t = varoom::funds::yield<int,std::string>(1);
    U v = varoom::funds::for_each<double>(t, [=](int i) {
        return varoom::funds::yield<double,std::string>(3.14 * i);
    });
    auto u = varoom::funds::parse(v, s);
    BOOST_CHECK_EQUAL(length(u), 1);
    BOOST_CHECK_EQUAL(access(u, 0).first, 3.14);
    BOOST_CHECK_EQUAL(access(u, 0).second, s);
}

BOOST_AUTO_TEST_CASE( fail1 )
{
    using T = varoom::funds::parser<std::string,int>;

    std::string s = "wibble";
    T t = varoom::funds::fail<int,std::string>();
    auto u = parse(t, s);
    BOOST_CHECK_EQUAL(length(u), 0);
}

BOOST_AUTO_TEST_CASE( one1 )
{
    using T = varoom::funds::parser<std::string,char>;

    std::string s = "wibble";
    T t = varoom::funds::one<std::string>();
    auto u = parse(t, s);
    BOOST_CHECK_EQUAL(length(u), 1);
    BOOST_CHECK_EQUAL(access(u, 0).first, 'w');
    BOOST_CHECK_EQUAL(access(u, 0).second, "ibble");
}

BOOST_AUTO_TEST_CASE( with1 )
{
    using V = std::pair<int,std::string>;
    using T = varoom::funds::parser<std::string,int>;

    std::string s = "wibble";
    T t = [](std::string st) {
        return make({V(1, st), V(2, st), V(3, st), V(4, st)});
    };

    auto u = parse(t, s);
    BOOST_CHECK_EQUAL(length(u), 4);
    BOOST_CHECK_EQUAL(access(u, 0).first, 1);
    BOOST_CHECK_EQUAL(access(u, 1).first, 2);
    BOOST_CHECK_EQUAL(access(u, 2).first, 3);
    BOOST_CHECK_EQUAL(access(u, 3).first, 4);

    T w = with(t, [](int x) {
        return bool((x & 1) == 0);
    });
    auto v = parse(w, s);
    BOOST_CHECK_EQUAL(length(v), 2);
    BOOST_CHECK_EQUAL(access(v, 0).first, 2);
    BOOST_CHECK_EQUAL(access(v, 1).first, 4);
}


BOOST_AUTO_TEST_CASE( orelse1 )
{
    using T = varoom::funds::parser<std::string,char>;

    std::string s = "wibble";
    T t = varoom::funds::one<std::string>();
    T u = with(t, [](char c) { return oneof("aeiou", c); });
    T v = with(t, [](char c) { return !oneof("aeiou", c); });
    T w = orelse(u, v);
    auto x = parse(w, s);
    BOOST_CHECK_EQUAL(length(x), 1);
}

