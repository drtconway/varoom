#include "varoom/util/subtext.hpp"

#define BOOST_TEST_MODULE subtext tests
#include <boost/test/unit_test.hpp>

BOOST_TEST_CASE( cons )
{
    using namespace std;
    using namespace varoom;

    string s("hello");
    subtext t1(s);
    BOOST_CHECK_EQUAL(t1.size(), s.size());

    subtext t2(s.begin(), s.end());
    BOOST_CHECK_EQUAL(t2.size(), s.size());

    string u = t2;
    BOOST_CHECK_EQUAL(u, s);
}

BOOST_TEST_CASE( split )
{
    using namespace std;
    using namespace varoom;

    string s("the quick  brown fox ");
    subtext t1(s);
    
    vector<subtext> parts1;
    t1.split(' ', parts);

    BOOST_CHECK_EQUAL(parts1.size(), 6);
    BOOST_CHECK_EQUAL(parts1[0], string("the"));
    BOOST_CHECK_EQUAL(parts1[1], string("quick"));
    BOOST_CHECK_EQUAL(parts1[2], string(""));
    BOOST_CHECK_EQUAL(parts1[3], string("brown"));
    BOOST_CHECK_EQUAL(parts1[4], string("fox"));
    BOOST_CHECK_EQUAL(parts1[5], string(""));
}
