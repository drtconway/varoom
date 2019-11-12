#include "varoom/util/lazy.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE lazy tests
#include <boost/test/unit_test.hpp>

using namespace std;
using namespace varoom;

namespace // anonymous
{
    int baz;

    string foo()
    {
        baz = 1;
        return "wibble";
    }

    class qux
    {
    public:
        qux()
        {
        }

        string quux() const
        {
            baz = 2;
            return "wombat";
        }
    };
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( simple )
{
    baz = 0;
    lazy<string> x(foo);
    BOOST_CHECK_EQUAL(baz, 0);
    BOOST_CHECK_EQUAL(x.get(), string("wibble"));
    BOOST_CHECK_EQUAL(baz, 1);
}

BOOST_AUTO_TEST_CASE( classy )
{
    baz = 0;
    qux q;
    lazy<string> x(bind(&qux::quux, q));
    BOOST_CHECK_EQUAL(baz, 0);
    BOOST_CHECK_EQUAL(x.get(), string("wombat"));
    BOOST_CHECK_EQUAL(baz, 2);
}

BOOST_AUTO_TEST_CASE( lammy )
{
    baz = 0;
    lazy<string> x([]() {
        baz = 3;
        return string("wangle");
    });
    BOOST_CHECK_EQUAL(baz, 0);
    BOOST_CHECK_EQUAL(x.get(), string("wangle"));
    BOOST_CHECK_EQUAL(baz, 3);
}

BOOST_AUTO_TEST_CASE( classy_lammy )
{
    baz = 0;
    qux q;
    lazy<string> x([q]() {
        return q.quux();
    });
    BOOST_CHECK_EQUAL(baz, 0);
    BOOST_CHECK_EQUAL(x.get(), string("wombat"));
    BOOST_CHECK_EQUAL(baz, 2);
}
