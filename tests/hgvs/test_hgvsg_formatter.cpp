#include "varoom/hgvs/hgvsg_formatter.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hgvs_positions tests
#include <boost/test/unit_test.hpp>

using namespace std;
using namespace varoom::hgvs;

BOOST_AUTO_TEST_CASE( testSub )
{
    string res;
    hgvsg_formatter fmt([&](const string& s) { res = s; });

    fmt.sub("chr13", 32972301, "G", "A");
    BOOST_CHECK_EQUAL(res, "chr13:g.32972301G>A");
}

BOOST_AUTO_TEST_CASE( testIns )
{
    string res;
    hgvsg_formatter fmt([&](const string& s) { res = s; });

    fmt.ins("chr13", 32972449, 32972450, "G");
    BOOST_CHECK_EQUAL(res, "chr13:g.32972449_32972450insG");
}

BOOST_AUTO_TEST_CASE( testDel1 )
{
    string res;
    hgvsg_formatter fmt([&](const string& s) { res = s; });

    fmt.del("chr13", 32972454, 32972454);
    BOOST_CHECK_EQUAL(res, "chr13:g.32972454del");
}

BOOST_AUTO_TEST_CASE( testDel2 )
{
    string res;
    hgvsg_formatter fmt([&](const string& s) { res = s; });

    fmt.del("chr13", 32972454, 32972455);
    BOOST_CHECK_EQUAL(res, "chr13:g.32972454_32972455del");
}

BOOST_AUTO_TEST_CASE( testDelIns1 )
{
    string res;
    hgvsg_formatter fmt([&](const string& s) { res = s; });

    fmt.delins("chr13", 32972453, 32972453, "AAAGG");
    BOOST_CHECK_EQUAL(res, "chr13:g.32972453delinsAAAGG");
}

BOOST_AUTO_TEST_CASE( testDelIns2 )
{
    string res;
    hgvsg_formatter fmt([&](const string& s) { res = s; });

    fmt.delins("chr13", 32972453, 32972457, "AAAGG");
    BOOST_CHECK_EQUAL(res, "chr13:g.32972453_32972457delinsAAAGG");
}

BOOST_AUTO_TEST_CASE( testDup1 )
{
    string res;
    hgvsg_formatter fmt([&](const string& s) { res = s; });

    fmt.dup("chr13", 32972454, 32972454);
    BOOST_CHECK_EQUAL(res, "chr13:g.32972454dup");
}

BOOST_AUTO_TEST_CASE( testDup2 )
{
    string res;
    hgvsg_formatter fmt([&](const string& s) { res = s; });

    fmt.dup("chr13", 32972454, 32972455);
    BOOST_CHECK_EQUAL(res, "chr13:g.32972454_32972455dup");
}

BOOST_AUTO_TEST_CASE( testInv )
{
    string res;
    hgvsg_formatter fmt([&](const string& s) { res = s; });

    fmt.inv("chr13", 32972454, 32972455);
    BOOST_CHECK_EQUAL(res, "chr13:g.32972454_32972455inv");
}

