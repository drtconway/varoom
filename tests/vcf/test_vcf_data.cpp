#include "varoom/vcf/vcf_data.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE vcf_data tests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( testInt )
{
    using namespace std;
    using namespace varoom::vcf;

    const int64_t i = 12345;
    const int64_t j = 54321;

    vcf_data vInt1(i);
    BOOST_CHECK_EQUAL(vInt1.type(), INT);
    BOOST_CHECK_EQUAL(vInt1.asInt(), i);

    vcf_data vInt2(vInt1);
    BOOST_CHECK_EQUAL(vInt2.type(), INT);
    BOOST_CHECK_EQUAL(vInt2.asInt(), i);

    vcf_data vInt3(j);
    BOOST_CHECK_EQUAL(vInt3.type(), INT);
    BOOST_CHECK_EQUAL(vInt3.asInt(), j);

    vInt3 = vInt2;
    BOOST_CHECK_EQUAL(vInt3.type(), INT);
    BOOST_CHECK_EQUAL(vInt3.asInt(), i);
}

