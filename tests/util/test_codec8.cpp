#include "varoom/util/codec8.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE codec8 tests
#include <boost/test/unit_test.hpp>

using namespace std;
using namespace varoom;
using namespace varoom::util;

namespace // anonymous
{
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( enc1 )
{
    uint64_t x = 0;
    vector<uint8_t> ys;
    size_t n = codec8::encode(x, ys);
    BOOST_CHECK_EQUAL(n, 1);
    BOOST_CHECK_EQUAL(ys.size(), 1);
    BOOST_CHECK_EQUAL(ys[0], 0);
}

BOOST_AUTO_TEST_CASE( enc2 )
{
    uint64_t x = 1;
    vector<uint8_t> ys;
    size_t n = codec8::encode(x, ys);
    BOOST_CHECK_EQUAL(n, 1);
    BOOST_CHECK_EQUAL(ys.size(), 1);
    BOOST_CHECK_EQUAL(ys[0], 2);
}

BOOST_AUTO_TEST_CASE( enc3 )
{
    uint64_t x = 123;
    vector<uint8_t> ys;
    size_t n = codec8::encode(x, ys);
    BOOST_CHECK_EQUAL(n, 1);
    BOOST_CHECK_EQUAL(ys.size(), 1);
    BOOST_CHECK_EQUAL(ys[0], 246);
}

BOOST_AUTO_TEST_CASE( enc4 )
{
    uint64_t x = 1234;
    vector<uint8_t> ys;
    size_t n = codec8::encode(x, ys);
    BOOST_CHECK_EQUAL(n, 2);
    BOOST_CHECK_EQUAL(ys.size(), 2);
    BOOST_CHECK_EQUAL(ys[0], 19);
    BOOST_CHECK_EQUAL(ys[1], 164);
}

BOOST_AUTO_TEST_CASE( dec1 )
{
    vector<uint8_t> ys{0};
    auto itr = ys.begin();
    uint64_t x = codec8::decode(itr);
    BOOST_CHECK_EQUAL(itr == ys.end(), true);
    BOOST_CHECK_EQUAL(x, 0);
}

BOOST_AUTO_TEST_CASE( dec2 )
{
    vector<uint8_t> ys{19, 164};
    auto itr = ys.begin();
    uint64_t x = codec8::decode(itr);
    BOOST_CHECK_EQUAL(itr == ys.end(), true);
    BOOST_CHECK_EQUAL(x, 1234);
}

