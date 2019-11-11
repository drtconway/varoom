#include "varoom/util/strong_typedef.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE strong_typedef tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
    struct meter : varoom::strong_typedef<meter, int>
    {
        using strong_typedef::strong_typedef;
    };

    struct inch : varoom::strong_typedef<inch, int>, varoom::integer_arithmetic<inch>
    {
        using strong_typedef::strong_typedef;
    };
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( test1 )
{
    meter m(3);
}

BOOST_AUTO_TEST_CASE( test2 )
{
    inch i(3);
    inch j(2);

    BOOST_CHECK_EQUAL(static_cast<int>(i + j), 5);
    BOOST_CHECK_EQUAL(static_cast<int>(i - j), 1);
    BOOST_CHECK_EQUAL(static_cast<int>(i * j), 6);
    BOOST_CHECK_EQUAL(static_cast<int>(i / j), 1);
    BOOST_CHECK_EQUAL(static_cast<int>(i % j), 1);
}
