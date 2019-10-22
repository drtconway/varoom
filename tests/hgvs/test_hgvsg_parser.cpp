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

BOOST_AUTO_TEST_CASE( testDel )
{
    using namespace std;
    using namespace varoom::hgvs;

    std::string res;
    hgvsg_formatter fmt([&](auto str){ res = str; });

    const std::string txt1 = "chr1:g.12345del";
    hgvsg_parser::parse(txt1, fmt);
    BOOST_CHECK_EQUAL(res, txt1);

    const std::string txt2 = "chr1:g.12345_12349del";
    hgvsg_parser::parse(txt2, fmt);
    BOOST_CHECK_EQUAL(res, txt2);

}

BOOST_AUTO_TEST_CASE( testDelIns )
{
    using namespace std;
    using namespace varoom::hgvs;

    std::string res;
    hgvsg_formatter fmt([&](auto str){ res = str; });

    const std::string txt1 = "chr1:g.12345delinsACGT";
    hgvsg_parser::parse(txt1, fmt);
    BOOST_CHECK_EQUAL(res, txt1);

    const std::string txt2 = "chr1:g.12345_12349delinsACGT";
    hgvsg_parser::parse(txt2, fmt);
    BOOST_CHECK_EQUAL(res, txt2);

}

BOOST_AUTO_TEST_CASE( testDup )
{
    using namespace std;
    using namespace varoom::hgvs;

    std::string res;
    hgvsg_formatter fmt([&](auto str){ res = str; });

    const std::string txt1 = "chr1:g.12345dup";
    hgvsg_parser::parse(txt1, fmt);
    BOOST_CHECK_EQUAL(res, txt1);

    const std::string txt2 = "chr1:g.12345_12349dup";
    hgvsg_parser::parse(txt2, fmt);
    BOOST_CHECK_EQUAL(res, txt2);

}

BOOST_AUTO_TEST_CASE( testInv )
{
    using namespace std;
    using namespace varoom::hgvs;

    std::string res;
    hgvsg_formatter fmt([&](auto str){ res = str; });

    const std::string txt = "chr1:g.12345_12346inv";
    hgvsg_parser::parse(txt, fmt);
    BOOST_CHECK_EQUAL(res, txt);
}

