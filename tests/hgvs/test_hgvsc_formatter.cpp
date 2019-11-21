#include "varoom/hgvs/hgvsc_formatter.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hgvsc_formatter tests
#include <boost/test/unit_test.hpp>

using namespace std;
using namespace varoom::hgvs;

namespace // anonymous
{
    hgvsc_locus l01(UTR5, hgvsc_position(-451), hgvsc_relative_position(0));
    hgvsc_locus l02(UTR5, hgvsc_position(-450), hgvsc_relative_position(0));
    hgvsc_locus l11(CODING, hgvsc_position(450), hgvsc_relative_position(0));
    hgvsc_locus l12(CODING, hgvsc_position(451), hgvsc_relative_position(0));
    hgvsc_locus l21(INTRON, hgvsc_position(450), hgvsc_relative_position(10));
    hgvsc_locus l22(INTRON, hgvsc_position(450), hgvsc_relative_position(11));
    hgvsc_locus l31(INTRON, hgvsc_position(450), hgvsc_relative_position(-11));
    hgvsc_locus l32(INTRON, hgvsc_position(450), hgvsc_relative_position(-10));
    hgvsc_locus l41(UTR3, hgvsc_position(450), hgvsc_relative_position(0));
    hgvsc_locus l42(UTR3, hgvsc_position(451), hgvsc_relative_position(0));
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( testSub )
{
    string res;
    hgvsc_formatter fmt([&](const string& s) { res = s; });

    {
        fmt.sub("NM_000059.3", l01, "G", "A");
        BOOST_CHECK_EQUAL(res, "NM_000059.3:c.-451G>A");
    }
    {
        fmt.sub("NM_000059.3", l11, "G", "A");
        BOOST_CHECK_EQUAL(res, "NM_000059.3:c.450G>A");
    }
    {
        fmt.sub("NM_000059.3", l21, "G", "A");
        BOOST_CHECK_EQUAL(res, "NM_000059.3:c.450+10G>A");
    }
    {
        fmt.sub("NM_000059.3", l31, "G", "A");
        BOOST_CHECK_EQUAL(res, "NM_000059.3:c.450-11G>A");
    }
    {
        fmt.sub("NM_000059.3", l41, "G", "A");
        BOOST_CHECK_EQUAL(res, "NM_000059.3:c.*450G>A");
    }
}

BOOST_AUTO_TEST_CASE( testIns )
{
    string res;
    hgvsc_formatter fmt([&](const string& s) { res = s; });

    {
        fmt.ins("NM_000059.3", l01, l02, "G");
        BOOST_CHECK_EQUAL(res, "NM_000059.3:c.-451_-450insG");
    }
}

BOOST_AUTO_TEST_CASE( testDel1 )
{
    string res;
    hgvsc_formatter fmt([&](const string& s) { res = s; });

    {
        fmt.del("NM_000059.3", l01, l01);
        BOOST_CHECK_EQUAL(res, "NM_000059.3:c.-451del");
    }
}

BOOST_AUTO_TEST_CASE( testDel2 )
{
    string res;
    hgvsc_formatter fmt([&](const string& s) { res = s; });

    {
        fmt.del("NM_000059.3", l01, l02);
        BOOST_CHECK_EQUAL(res, "NM_000059.3:c.-451_-450del");
    }
}

BOOST_AUTO_TEST_CASE( testDelIns1 )
{
    string res;
    hgvsc_formatter fmt([&](const string& s) { res = s; });

    {
        fmt.delins("NM_000059.3", l01, l01, "G");
        BOOST_CHECK_EQUAL(res, "NM_000059.3:c.-451delinsG");
    }
}

BOOST_AUTO_TEST_CASE( testDelIns2 )
{
    string res;
    hgvsc_formatter fmt([&](const string& s) { res = s; });

    {
        fmt.delins("NM_000059.3", l01, l02, "G");
        BOOST_CHECK_EQUAL(res, "NM_000059.3:c.-451_-450delinsG");
    }
}

BOOST_AUTO_TEST_CASE( testDup1 )
{
    string res;
    hgvsc_formatter fmt([&](const string& s) { res = s; });

    {
        fmt.dup("NM_000059.3", l01, l01);
        BOOST_CHECK_EQUAL(res, "NM_000059.3:c.-451dup");
    }
}

BOOST_AUTO_TEST_CASE( testDup2 )
{
    string res;
    hgvsc_formatter fmt([&](const string& s) { res = s; });

    {
        fmt.dup("NM_000059.3", l01, l02);
        BOOST_CHECK_EQUAL(res, "NM_000059.3:c.-451_-450dup");
    }
}

BOOST_AUTO_TEST_CASE( testInv )
{
    string res;
    hgvsc_formatter fmt([&](const string& s) { res = s; });

    {
        fmt.inv("NM_000059.3", l01, l02);
        BOOST_CHECK_EQUAL(res, "NM_000059.3:c.-451_-450inv");
    }
}

