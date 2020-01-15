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

    template<typename T, typename U>
    parser<std::pair<T,U>> cons(parser<T> p, parser<U> q)
    {
        return (p >>= [=](T t) {
            return (q >>= [=](U u) {
                return yield<parser,std::pair<T,U>>(std::make_pair(t, u));
            });
        });
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
