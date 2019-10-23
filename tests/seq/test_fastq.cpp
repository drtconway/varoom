#include "varoom/seq/fastq.hpp"

#include <initializer_list>
#include <sstream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE fastq tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
    std::string compose(std::initializer_list<const char*> p_parts)
    {
        std::string s;
        for (auto i = p_parts.begin(); i != p_parts.end(); ++i)
        {
            s.insert(s.size(), *i);
        }
        return s;
    }
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( test_parser )
{
    using namespace std;
    using namespace varoom::seq;

    const std::string fq = compose({
        "@NS500817:604:HLY3VBGXC:1:11101:23926:2283 1:N:0:TAGGATGA\n",
        "CCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGA\n",
        "+\n",
        "AA/A6EAEE6E66A/E6EE/////EEEE6AEE66AE/EEEEE//AEE6A///AAEAEAEE/////A/EEE/E/E/\n",
        "@NS500817:604:HLY3VBGXC:1:11101:4317:2283 1:N:0:TAGGATGA\n",
        "AGCTGGGCTGGCCTGAGCTGGGCTGGGCGAGGCTGGGCTGGGCTGGGCTGGGCTTGGACGAGCTGAGATCGGAAG\n",
        "+\n",
        "AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<E\n",
        "@NS500817:604:HLY3VBGXC:1:11101:10495:2283 1:N:0:TAGGATGA\n",
        "TGAGACAAGAGGAAGGCATCTGTCTCCTGCCCTTCCCTGGGCAATGGAATGTCTCGGTGTAAAACCCGATTGTAT\n",
        "+\n",
        "AAAAAEEEEEEEEEEEEEEEEEEEEAAEEEE/EEAAEEEEEEEEEEEEEEEAEEAEEEEEEE6EEE<EAEEEEEE\n",
        "@NS500817:604:HLY3VBGXC:1:11101:26104:2284 1:N:0:TAGGATGA\n",
        "GACCTGGATGGTTCTCCCCGTCGCGCTCTGCCCAGCTCTCCCTCCCTCCGTCCCTCCCTAGATCGGAAGAGCACA\n",
        "+\n",
        "A<AAAEEEEEEEEEEAEEEAEAAEEEEEEEEAEAEEEAEEEEEEEEEEEEEEEEAAEEEEAEEAEEEEEEEEEEE\n"
    });

    istringstream in(fq);
    fastq_reader r(in);
    BOOST_CHECK_EQUAL(r.more(), true);
    BOOST_CHECK_EQUAL(std::get<0>(*r), "NS500817:604:HLY3VBGXC:1:11101:23926:2283 1:N:0:TAGGATGA");
    BOOST_CHECK_EQUAL(std::get<1>(*r), "CCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGA");
    BOOST_CHECK_EQUAL(std::get<2>(*r), "");
    BOOST_CHECK_EQUAL(std::get<3>(*r), "AA/A6EAEE6E66A/E6EE/////EEEE6AEE66AE/EEEEE//AEE6A///AAEAEAEE/////A/EEE/E/E/");
    ++r;
    BOOST_CHECK_EQUAL(r.more(), true);
    BOOST_CHECK_EQUAL(std::get<0>(*r), "NS500817:604:HLY3VBGXC:1:11101:4317:2283 1:N:0:TAGGATGA");
    ++r;
    BOOST_CHECK_EQUAL(r.more(), true);
    BOOST_CHECK_EQUAL(std::get<0>(*r), "NS500817:604:HLY3VBGXC:1:11101:10495:2283 1:N:0:TAGGATGA");
    ++r;
    BOOST_CHECK_EQUAL(r.more(), true);
    BOOST_CHECK_EQUAL(std::get<0>(*r), "NS500817:604:HLY3VBGXC:1:11101:26104:2284 1:N:0:TAGGATGA");
    ++r;
    BOOST_CHECK_EQUAL(r.more(), false);
}

