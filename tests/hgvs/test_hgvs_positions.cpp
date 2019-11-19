#include "varoom/hgvs/hgvs_positions.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hgvs_positions tests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( testGenomic )
{
    using namespace std;
    using namespace varoom::hgvs;

    genomic_locus g0(12345);
    BOOST_CHECK_EQUAL(static_cast<uint64_t>(g0), 12345);

    hgvsg_locus g1(12345);
    BOOST_CHECK_EQUAL(static_cast<uint64_t>(g1), 12345);

    tx_position g2(12345);
    BOOST_CHECK_EQUAL(static_cast<int64_t>(g2), 12345);

    hgvsc_position g3(12345);
    BOOST_CHECK_EQUAL(static_cast<int64_t>(g3), 12345);
}
