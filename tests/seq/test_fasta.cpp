#include "varoom/seq/fasta.hpp"

#include <initializer_list>
#include <sstream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE fasta tests
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

    const std::string fa = compose({
        ">foo\n",
        "AGCCAGTTCCGGGTGTCGCCGCTGGATCGGACCTGGAACCTGGGCGAGACAGTGGAGCTG\n",
        ">bar\n",
        "AAGTGCCAGGTGCTGCTGTCCAACCCGACGTCGGGCTGCTCGTGGCTCTTCCAGCCGCGC\n",
        "GGCGCCGCCGCCAGTCCCACCTTCCTCCTATACCTCTCCCAAAACAAGCCCAAGGCGGCC\n",
        "GAGGGGCTGGACACCCAGCGGTTCTCGGGCAAGAGGTTGGGGGACACCTTCGTCCTCACC\n",
        ">baz\n",
        ">qux\n",
        "CTGAGCGACTTCCGCCGAGAGAACGAGGGCTACTATTTCTGCTCGGCCCTGAGCAACTCC\n",
        "ATCATGTACTTCAGCCACTTCGTGCCGGTCTTCCTGCCA\n"
    });

    istringstream in(fa);
    fasta_reader r(in);
    BOOST_CHECK_EQUAL(r.more(), true);
    BOOST_CHECK_EQUAL((*r).first, "foo");
    BOOST_CHECK_EQUAL((*r).second, "AGCCAGTTCCGGGTGTCGCCGCTGGATCGGACCTGGAACCTGGGCGAGACAGTGGAGCTG");
    ++r;
    BOOST_CHECK_EQUAL(r.more(), true);
    BOOST_CHECK_EQUAL((*r).first, "bar");
    BOOST_CHECK_EQUAL((*r).second.size(), 180);
    ++r;
    BOOST_CHECK_EQUAL(r.more(), true);
    BOOST_CHECK_EQUAL((*r).first, "baz");
    BOOST_CHECK_EQUAL((*r).second.size(), 0);
    ++r;
    BOOST_CHECK_EQUAL(r.more(), true);
    BOOST_CHECK_EQUAL((*r).first, "qux");
    BOOST_CHECK_EQUAL((*r).second.size(), 99);
    ++r;
    BOOST_CHECK_EQUAL(r.more(), false);
}

