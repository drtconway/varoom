#include "varoom/hgvs/hgvsc_locus.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hgvsc_locus tests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( testUtr5 )
{
    using namespace std;
    using namespace varoom::hgvs;

    {
        hgvsc_locus l(UTR5, hgvsc_position(-5), hgvsc_relative_position(0));
        BOOST_CHECK_EQUAL(l.str(), "-5");
    }
    {
        hgvsc_locus l(UTR5, hgvsc_position(-5), hgvsc_relative_position(-7));
        BOOST_CHECK_EQUAL(l.str(), "-5-7");
    }
    {
        hgvsc_locus l(UTR5, hgvsc_position(-6), hgvsc_relative_position(7));
        BOOST_CHECK_EQUAL(l.str(), "-6+7");
    }
}

BOOST_AUTO_TEST_CASE( testUtr3 )
{
    using namespace std;
    using namespace varoom::hgvs;

    {
        hgvsc_locus l(UTR3, hgvsc_position(5), hgvsc_relative_position(0));
        BOOST_CHECK_EQUAL(l.str(), "*5");
    }
    {
        hgvsc_locus l(UTR3, hgvsc_position(5), hgvsc_relative_position(7));
        BOOST_CHECK_EQUAL(l.str(), "*5+7");
    }
    {
        hgvsc_locus l(UTR3, hgvsc_position(6), hgvsc_relative_position(-7));
        BOOST_CHECK_EQUAL(l.str(), "*6-7");
    }
}

BOOST_AUTO_TEST_CASE( testInterAlia )
{
    using namespace std;
    using namespace varoom::hgvs;

    {
        hgvsc_locus l(CODING, hgvsc_position(5), hgvsc_relative_position(0));
        BOOST_CHECK_EQUAL(l.str(), "5");
    }
    {
        hgvsc_locus l(INTRON, hgvsc_position(5), hgvsc_relative_position(7));
        BOOST_CHECK_EQUAL(l.str(), "5+7");
    }
    {
        hgvsc_locus l(INTRON, hgvsc_position(6), hgvsc_relative_position(-7));
        BOOST_CHECK_EQUAL(l.str(), "6-7");
    }
}
