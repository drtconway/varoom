#include "varoom/hgvs/hgvs_positions.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hgvs_positions tests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( testGenomic )
{
    using namespace std;
    using namespace varoom::hgvs;

    genomic_locus g0(12345);
    BOOST_CHECK_EQUAL(g0(), 12345);
    hgvsg_locus g1 = static_cast<hgvsg_locus>(g0);
    BOOST_CHECK_EQUAL(g1(), 12346);
    genomic_locus g2 = static_cast<genomic_locus>(g1);
    BOOST_CHECK_EQUAL(g2(), 12345);

}
