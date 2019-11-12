#include "varoom/vcf/vcf_writer.hpp"
#include "varoom/vcf/vcf_parser.hpp"

#include <fstream>
#include <sstream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE vcf_writer tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( test1 )
{
    using namespace std;
    using namespace varoom::vcf;

    ifstream in("tests/data/exac-few.vcf");
    ostringstream out;

    vcf_writer T(out);
    vcf_parser P(T);
    P.parse(in);

    string s = out.str();
}
