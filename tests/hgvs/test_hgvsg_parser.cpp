#include "varoom/hgvs/hgvsg_parser.hpp"
#include "varoom/hgvs/hgvsg_formatter.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hgvsg_parser tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
    class string_capture
    {
    public:
        std::string str;

        void accept(const std::string& p_str)
        {
            str = p_str;
        }
    };
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( testSub )
{
    using namespace std;
    using namespace varoom::hgvs;

    std::string res;
    hgvsg_formatter fmt([&](auto str){ res = str; });

    const std::string txt = "chr1:g.12345A>C";
    hgvsg_parser::parse(txt, fmt);
    BOOST_CHECK_EQUAL(res, txt);
}

BOOST_AUTO_TEST_CASE( testIns )
{
    using namespace std;
    using namespace varoom::hgvs;

    std::string res;
    hgvsg_formatter fmt([&](auto str){ res = str; });

    const std::string txt = "chr1:g.12345_12346insACGT";
    hgvsg_parser::parse(txt, fmt);
    BOOST_CHECK_EQUAL(res, txt);
}

