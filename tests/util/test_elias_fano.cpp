#include "varoom/util/elias_fano.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ellias_fano tests
#include <boost/test/unit_test.hpp>

#include <random>
#include <iostream>

using namespace std;
using namespace varoom;

namespace // anonymous
{

}
// namespace anonymous

BOOST_AUTO_TEST_CASE( enc0 )
{
    std::vector<uint64_t> X({(1ULL << 16) + 3, 1361038616ULL});
    std::vector<uint64_t> N;
    std::vector<uint64_t> Y;
    for (size_t i = 0; i < X.size(); ++i)
    {
        uint64_t x = X[i];
        size_t n = varoom::elias::ilog2c(x);
        N.push_back(n);
        uint64_t y = varoom::elias::rev(n, x);
        BOOST_REQUIRE_EQUAL(y & 1, 1);
        Y.push_back(y);
    }

    bitstream B;
    for (uint64_t i = 0; i < X.size(); ++i)
    {
        uint64_t x = X[i];
        varoom::elias::encode(B, x);
    }

    varoom::bitstream::iterator itr = B.begin();
    varoom::bitstream::iterator jtr = B.begin();
    {
        size_t t = 0;
        uint64_t x;

        x = B.at(N[0] - 1, t);
        BOOST_REQUIRE_EQUAL(x, 0);
        BOOST_REQUIRE_EQUAL(itr.pop_front(N[0] - 1), 0);
        BOOST_REQUIRE_EQUAL(varoom::elias::count_leading_zeros(jtr), N[0] - 1);
        t += N[0] - 1;

        BOOST_CHECK_EQUAL(t, 16);
        BOOST_CHECK_EQUAL(itr.pos, 16);
        BOOST_CHECK_EQUAL(jtr.pos, 16);

        x = B.at(N[0], t);
        BOOST_REQUIRE_EQUAL(x, Y[0]);
        BOOST_REQUIRE_EQUAL(itr.pop_front(N[0]), Y[0]);
        BOOST_REQUIRE_EQUAL(jtr.pop_front(N[0]), Y[0]);
        t += N[0];

        BOOST_CHECK_EQUAL(t, 33);
        BOOST_CHECK_EQUAL(itr.pos, 33);
        BOOST_CHECK_EQUAL(jtr.pos, 33);

        x = B.at(N[1] - 1, t);
        BOOST_REQUIRE_EQUAL(x, 0);
        BOOST_REQUIRE_EQUAL(itr.pop_front(N[1] - 1), 0);
        BOOST_REQUIRE_EQUAL(varoom::elias::count_leading_zeros(jtr), N[1] - 1);
        t += N[1] - 1;

        BOOST_CHECK_EQUAL(t, 63);
        BOOST_CHECK_EQUAL(itr.pos, 63);
        BOOST_CHECK_EQUAL(jtr.pos, 63);

        x = B.at(N[1], t);
        BOOST_REQUIRE_EQUAL(x, Y[1]);
        BOOST_REQUIRE_EQUAL(itr.pop_front(N[1]), Y[1]);
        BOOST_REQUIRE_EQUAL(jtr.pop_front(N[1]), Y[1]);
        t += N[1];
    }

    {
        itr = B.begin();
        for (size_t i = 0; i < B.size(); ++i, ++itr)
        {
            BOOST_REQUIRE_EQUAL(*itr, B.at(1, i));
        }
    }
    
    {
        itr = B.begin();
        for (size_t i = 0; i < B.size(); ++i)
        {
            BOOST_REQUIRE_EQUAL(itr.pop_front(1), B.at(1, i));
        }
    }
    
    itr = B.begin();
    uint64_t x;

    x = varoom::elias::decode(itr);
    BOOST_REQUIRE_EQUAL(x, X[0]);

    BOOST_CHECK_EQUAL(itr.pos, 33);

    size_t n = varoom::elias::count_leading_zeros(itr);
    BOOST_REQUIRE_EQUAL(n, N[1] - 1);
}

BOOST_AUTO_TEST_CASE( enc1 )
{
    for (uint64_t x = 1; x < 1000; ++x)
    {
        bitstream B;
        varoom::elias::encode(B, x);
        bitstream::iterator itr = B.begin();
        uint64_t y = varoom::elias::decode(itr);
        BOOST_CHECK_EQUAL(y, x);
    }
}

BOOST_AUTO_TEST_CASE( enc2 )
{
    std::vector<uint64_t> X;
    bitstream B;

    std::mt19937_64 rng(17);
    std::uniform_int_distribution<uint64_t> Z(1, (1ULL << 32) - 1);
    for (uint64_t i = 0; i < 1024*1024; ++i)
    {
        uint64_t x = Z(rng);
        X.push_back(x);
        varoom::elias::encode(B, x);
    }

    varoom::bitstream::iterator itr = B.begin();
    size_t t = 0;
    for (uint64_t i = 0; i < X.size(); ++i)
    {
        BOOST_REQUIRE_EQUAL(itr.pos, t);
        uint64_t x = X[i];
        uint64_t y = varoom::elias::decode(itr);
        if (x != y)
        {
            std::cerr << i
                << '\t' << itr.pos
                << '\t' << (itr.pos % 64)
                << '\t' << x
                << '\t' << varoom::elias::ilog2c(x)
                << '\t' << y
                << std::endl;
        }
        BOOST_REQUIRE_EQUAL(y, x);
        size_t n = varoom::elias::ilog2c(x);
        t += 2*n - 1;
    }
}
