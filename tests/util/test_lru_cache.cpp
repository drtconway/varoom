#include "varoom/util/lru_cache.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE lru_cache tests
#include <boost/test/unit_test.hpp>

using namespace std;
using namespace varoom;

namespace // anonymous
{
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( simple_1 )
{
    int got = 0;

    varoom::lru_cache<int,int> C(10, [&](const int& p_key) { got = p_key; return -2*p_key; });

    for (int i = 0; i < 10; ++i)
    {
        int r = C[i];
        BOOST_CHECK_EQUAL(got, i);
        BOOST_CHECK_EQUAL(r, -2*i);
    }

}
