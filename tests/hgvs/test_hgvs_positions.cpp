#include "varoom/hgvs/hgvs_positions.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hgvs_positions tests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( testGenomic )
{
    using namespace std;
    using namespace varoom::hgvs;

    gpos_0 g0(12345);
    BOOST_CHECK_EQUAL(g0(), 12345);
    gpos_1 g1 = static_cast<gpos_1>(g0);
    BOOST_CHECK_EQUAL(g1(), 12346);
    gpos_0 g2 = static_cast<gpos_0>(g1);
    BOOST_CHECK_EQUAL(g2(), 12345);

}
