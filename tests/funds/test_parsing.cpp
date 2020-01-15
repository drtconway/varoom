#include "varoom/funds/parsing.hpp"

#include <iostream>

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

    template<typename T, typename U>
    parser<std::pair<T,U>> cons(parser<T> p, parser<U> q)
    {
        return (p >>= [=](T t) {
            return (q >>= [=](U u) {
                return yield<parser,std::pair<T,U>>(std::make_pair(t, u));
            });
        });
    }

    bool vowel(char c)
    {
        switch (c)
        {
            case 'A':
            case 'a':
            case 'E':
            case 'e':
            case 'I':
            case 'i':
            case 'O':
            case 'o':
            case 'U':
            case 'u':
            {
                return true;
            }
            default:
            {
                return false;
            }
        }
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

BOOST_AUTO_TEST_CASE( sym0 )
{
    using P = parser<char>;
    P p = sym();
    auto r = parse(p, symbols(""));
    BOOST_CHECK_EQUAL(length(r), 0);
}

BOOST_AUTO_TEST_CASE( sym1 )
{
    using P = parser<char>;
    P p = sym();
    auto r = parse(p, symbols("a"));
    BOOST_CHECK_EQUAL(length(r), 1);
    BOOST_CHECK_EQUAL(access(r, 0).first, 'a');
    symbols t = access(r, 0).second;
    std::string u(t.begin(), t.end());
    BOOST_CHECK_EQUAL(u, "");
}

BOOST_AUTO_TEST_CASE( sym2a )
{
    using P = parser<char>;
    using Q = parser<std::pair<char,char>>;
    P p = sym();
    Q q = cons(p, p);
    auto r = parse(q, symbols("a"));
    BOOST_CHECK_EQUAL(length(r), 0);
}

BOOST_AUTO_TEST_CASE( sym2b )
{
    using P = parser<char>;
    using Q = parser<std::pair<char,char>>;
    P p = sym();
    Q q = cons(p, p);
    auto r = parse(q, symbols("ab"));
    BOOST_CHECK_EQUAL(length(r), 1);
    BOOST_CHECK_EQUAL(access(r, 0).first.first, 'a');
    BOOST_CHECK_EQUAL(access(r, 0).first.second, 'b');
    symbols t = access(r, 0).second;
    std::string u(t.begin(), t.end());
    BOOST_CHECK_EQUAL(u, "");
}

BOOST_AUTO_TEST_CASE( sat1 )
{
    using P = parser<char>;
    P p = sym();
    P q = sat(p, [](char c) { return vowel(c); });
    auto r = parse(q, symbols("ab"));
    BOOST_CHECK_EQUAL(length(r), 1);
    BOOST_CHECK_EQUAL(access(r, 0).first, 'a');
    symbols t = access(r, 0).second;
    std::string u(t.begin(), t.end());
    BOOST_CHECK_EQUAL(u, "b");
}

BOOST_AUTO_TEST_CASE( sat2 )
{
    using P = parser<char>;
    P p = sym();
    P q = sat(p, [](char c) { return vowel(c); });
    auto r = parse(q, symbols("ba"));
    BOOST_CHECK_EQUAL(length(r), 0);
}

BOOST_AUTO_TEST_CASE( or1 )
{
    using P = parser<char>;
    P p = sym('a');
    P q = sym('b');
    P r = append(p, q);
    {
        auto s = parse(r, symbols(""));
        BOOST_CHECK_EQUAL(length(s), 0);
    }
    {
        auto s = parse(r, symbols("a"));
        BOOST_REQUIRE_EQUAL(length(s), 1);
        BOOST_CHECK_EQUAL(access(s, 0).first, 'a');
    }
    {
        auto s = parse(r, symbols("b"));
        BOOST_REQUIRE_EQUAL(length(s), 1);
        BOOST_CHECK_EQUAL(access(s, 0).first, 'b');
    }
}

BOOST_AUTO_TEST_CASE( many0_1 )
{
    using P = parser<char>;
    using Q = parser<list<char>>;
    P p = sym();
    Q q = many0(p);
    {
        auto s = parse(q, symbols(""));
        BOOST_REQUIRE_EQUAL(length(s), 1);
        BOOST_CHECK_EQUAL(access(s, 0).first.empty(), true);
    }
    {
        auto s = parse(q, symbols("a"));
        BOOST_REQUIRE_EQUAL(length(s), 1);
        BOOST_CHECK_EQUAL(length(access(s, 0).first), 1);
    }
    {
        auto s = parse(q, symbols("aa"));
        BOOST_REQUIRE_EQUAL(length(s), 1);
        BOOST_CHECK_EQUAL(length(access(s, 0).first), 2);
    }
    {
        auto s = parse(q, symbols("aaa"));
        BOOST_REQUIRE_EQUAL(length(s), 1);
        BOOST_CHECK_EQUAL(length(access(s, 0).first), 3);
    }
}

BOOST_AUTO_TEST_CASE( many1_1 )
{
    using P = parser<char>;
    using Q = parser<list<char>>;
    P p = sym();
    Q q = many1(p);
    {
        auto s = parse(q, symbols(""));
        BOOST_REQUIRE_EQUAL(length(s), 0);
    }
    {
        auto s = parse(q, symbols("a"));
        BOOST_REQUIRE_EQUAL(length(s), 1);
        BOOST_CHECK_EQUAL(length(access(s, 0).first), 1);
    }
    {
        auto s = parse(q, symbols("aa"));
        BOOST_REQUIRE_EQUAL(length(s), 1);
        BOOST_CHECK_EQUAL(length(access(s, 0).first), 2);
    }
    {
        auto s = parse(q, symbols("aaa"));
        BOOST_REQUIRE_EQUAL(length(s), 1);
        BOOST_CHECK_EQUAL(length(access(s, 0).first), 3);
    }
}

BOOST_AUTO_TEST_CASE( fmap1 )
{
    using P = parser<list<char>>;
    using Q = parser<std::string>;
    P p = many0(sym());
    static_assert(std::is_convertible<decltype(list_to_string), std::function<std::string(list<char>)>>::value);
    Q q = fmap(list_to_string, p);
    {
        auto s = parse(q, symbols("aaa"));
        BOOST_REQUIRE_EQUAL(length(s), 1);
        BOOST_CHECK_EQUAL(access(s, 0).first.size(), 3);
        BOOST_CHECK_EQUAL(access(s, 0).first, "aaa");
    }
}
