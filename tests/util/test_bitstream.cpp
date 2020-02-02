#include "varoom/util/bitstream.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE bitstream tests
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <random>

using namespace std;
using namespace varoom;

namespace // anonymous
{
    std::string bin(size_t n, uint64_t x)
    {
        std::string s;
        for (size_t i = 0; i < n; ++i)
        {
            s.push_back("01"[x&1]);
            x >>= 1;
        }
        return s;
    }

}
// namespace anonymous

BOOST_AUTO_TEST_CASE( enc00 )
{
    std::string S;
    bitstream B;

    B.push_back(16, 0ULL);
    S += bin(16, 0ULL);

    BOOST_CHECK_EQUAL(B.size(), 16);
    BOOST_CHECK_EQUAL(B.str(), S);

    B.push_back(17, 98305ULL);
    S += bin(17, 98305ULL);

    BOOST_CHECK_EQUAL(B.size(), 33);
    BOOST_CHECK_EQUAL(B.str(), S);

    B.push_back(30, 0ULL);
    S += bin(30, 0ULL);

    BOOST_CHECK_EQUAL(B.size(), 63);
    BOOST_CHECK_EQUAL(B.str(), S);

    B.push_back(31, 207223877ULL);
    S += bin(31, 207223877ULL);

    BOOST_CHECK_EQUAL(B.size(), 94);
    BOOST_CHECK_EQUAL(B.str(), S);

    BOOST_CHECK_EQUAL(B.at(16, 0), 0ULL);
    BOOST_CHECK_EQUAL(B.at(17, 16), 98305ULL);
    BOOST_CHECK_EQUAL(B.at(30, 33), 0ULL);
    BOOST_CHECK_EQUAL(B.at(31, 63), 207223877ULL);
}

BOOST_AUTO_TEST_CASE( enc0 )
{
    varoom::bitstream B;
    B.push_back(55, 0);
    B.push_back(11, 10);
    auto itr = B.begin();
    BOOST_REQUIRE_EQUAL(itr.pop_front(55), 0);
    uint64_t x = itr.pop_front(11);
    BOOST_REQUIRE_EQUAL(x, 10);
}

BOOST_AUTO_TEST_CASE( enc1 )
{
    varoom::bitstream B;
    uint64_t t = 0;
    for (uint64_t i = 0; i < 1000; ++i)
    {
        uint64_t n = 1 + (i % 64);
        uint64_t m = (1ULL << n) - 1;
        uint64_t x = i & m;
        B.push_back(n, x);
        t += n;
        BOOST_REQUIRE_EQUAL(B.size(), t);
    }
    BOOST_REQUIRE_EQUAL(B.size(), t);
    auto itr = B.begin();
    for (uint64_t i = 0; i < 1000; ++i)
    {
        uint64_t n = 1 + (i % 64);
        uint64_t m = (1ULL << n) - 1;
        uint64_t x = i & m;
        uint64_t y = itr.pop_front(n);
        BOOST_REQUIRE_EQUAL(y, x);
    }
}

BOOST_AUTO_TEST_CASE( enc2 )
{
    std::mt19937_64 rng(17);

    using item = std::pair<size_t,uint64_t>;
    std::vector<item> V;
    varoom::bitstream B;

    size_t t = 0;
    for (uint64_t i = 0; i < 10000; ++i)
    {
        uint64_t n = 1 + (i % 64);
        uint64_t m = (1ULL << n) - 1;
        std::uniform_int_distribution<uint64_t> Z(0, m);
        uint64_t x = Z(rng);
        V.push_back(item(n, x));
        B.push_back(n, x);
        t += n;
        BOOST_REQUIRE_EQUAL(B.size(), t);
    }
    BOOST_REQUIRE_EQUAL(B.size(), t);
    auto itr = B.begin();
    for (uint64_t i = 0; i < V.size(); ++i)
    {
        uint64_t n = V[i].first;
        uint64_t x = V[i].second;
        uint64_t y = itr.pop_front(n);
        BOOST_REQUIRE_EQUAL(y, x);
    }
}
