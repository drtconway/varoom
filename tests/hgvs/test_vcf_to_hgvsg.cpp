#include "varoom/hgvs/vcf_to_hgvsg.hpp"
#include "varoom/hgvs/hgvsg_formatter.hpp"
#include "varoom/vcf/vcf_parser.hpp"

#include <fstream>
#include <sstream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE vcf_to_hgvsg tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( testFile )
{
    using namespace std;
    using namespace varoom::hgvs;
    using namespace varoom::vcf;

    ifstream in("tests/data/exac-few.vcf");
    
    vector<string> hgvsgs;

    hgvsg_formatter F([&](auto& p_str) { hgvsgs.push_back(p_str); });
    vcf_to_hgvsg V(F);
    vcf_parser P(V);
    P.parse(in);
    for (size_t i = 0; i < hgvsgs.size(); ++i)
    {
        //cout << i << '\t' << hgvsgs[i] << endl;
    }
}
