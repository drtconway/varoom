#include "varoom/util/rope.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE rope tests
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <random>

namespace // anonymous
{
    const size_t B = 10;

    std::mt19937_64 rng;

    std::map<size_t, std::string> X;

    std::string load_block(size_t p_n)
    {
        std::uniform_int_distribution<char> U('A', 'Z');
        std::string x;
        rng.seed(p_n);
        for (size_t i = 0; i < (1ULL << B); ++i)
        {
            x.push_back(U(rng));
        }
        X[p_n] = x;
        return x;
    }
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( cons1 )
{
    varoom::rope r("foo");
    BOOST_CHECK_EQUAL(r.size(), 3);
    BOOST_CHECK_EQUAL(r.str(), "foo");
}

BOOST_AUTO_TEST_CASE( cons2 )
{
    varoom::rope r1("foo");
    varoom::rope r2("baz");
    varoom::rope r(r1, r2);
    BOOST_CHECK_EQUAL(r.size(), 6);
    BOOST_CHECK_EQUAL(r.str(), "foobaz");
}

BOOST_AUTO_TEST_CASE( cons2b )
{
    varoom::rope r1("foo");
    varoom::rope r2("baz");
    varoom::rope r = r1 + r2;
    BOOST_CHECK_EQUAL(r.size(), 6);
    BOOST_CHECK_EQUAL(r.str(), "foobaz");
}

BOOST_AUTO_TEST_CASE( cons3 )
{
    varoom::rope r1("foobarbaz");
    varoom::rope r(r1, 3, 6);
    BOOST_CHECK_EQUAL(r.size(), 3);
    BOOST_CHECK_EQUAL(r.str(), "bar");
}

BOOST_AUTO_TEST_CASE( cons4 )
{
    varoom::rope r1("foobarbaz");
    varoom::rope r2("gppcbscb0");
    varoom::rope r3(r1, r2);
    {
        varoom::rope r(r3, 3, 6);
        BOOST_CHECK_EQUAL(r.size(), 3);
        BOOST_CHECK_EQUAL(r.str(), "bar");
    }
    {
        varoom::rope r(r3, 6, 12);
        BOOST_CHECK_EQUAL(r.size(), 6);
        BOOST_CHECK_EQUAL(r.str(), "bazgpp");
    }
    {
        varoom::rope r(r3, 12, 15);
        BOOST_CHECK_EQUAL(r.size(), 3);
        BOOST_CHECK_EQUAL(r.str(), "cbs");
    }
}

BOOST_AUTO_TEST_CASE( cons5 )
{
    const size_t N = 100000;
    std::function<std::string(size_t)> f = load_block;
    varoom::rope r1(N, B, f);

    varoom::rope r = r1.slice(2000, 4000);
    BOOST_CHECK_EQUAL(X.size(), 0);
    BOOST_CHECK_EQUAL(r.size(), 2000);
    std::string s = r.str();
    BOOST_CHECK_EQUAL(X.size(), 3);
}
