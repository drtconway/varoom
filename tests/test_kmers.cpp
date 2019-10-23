#include "varoom/kmers.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE fasta tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( test_make_1 )
{
    using namespace std;
    using namespace varoom::kmers;

    const string s = "AGCCAGTTCCGGGTGTCGCCGCTGGATCGGACCTGGAACCTGGGCGAGACAGTGGAGCTG";
    const size_t K = 5;
    vector<kmer> v;
    make(s, K, v);
    BOOST_CHECK_EQUAL(v.size(), 56);
    BOOST_CHECK_EQUAL(render(K, v[0]), "AGCCA");
    BOOST_CHECK_EQUAL(render(K, v[5]), "GTTCC");
    BOOST_CHECK_EQUAL(render(K, v[10]), "GGGTG");
    BOOST_CHECK_EQUAL(render(K, v[15]), "TCGCC");
    BOOST_CHECK_EQUAL(render(K, v[20]), "GCTGG");
    BOOST_CHECK_EQUAL(render(K, v[25]), "ATCGG");
    BOOST_CHECK_EQUAL(render(K, v[30]), "ACCTG");
    BOOST_CHECK_EQUAL(render(K, v[35]), "GAACC");
    BOOST_CHECK_EQUAL(render(K, v[40]), "TGGGC");
    BOOST_CHECK_EQUAL(render(K, v[45]), "GAGAC");
    BOOST_CHECK_EQUAL(render(K, v[50]), "AGTGG");
    BOOST_CHECK_EQUAL(render(K, v[55]), "AGCTG");
}

BOOST_AUTO_TEST_CASE( test_make_2 )
{
    using namespace std;
    using namespace varoom::kmers;

    const string s = "AGCCAGTTCCGGGTGTCGCCGCTGGATCGGNCCTGGAACCTGGGCGAGACAGTGGAGCTG";
    const size_t K = 5;
    vector<kmer> v;
    make(s, K, v);
    BOOST_CHECK_EQUAL(v.size(), 51);
    BOOST_CHECK_EQUAL(render(K, v[0]), "AGCCA");
    BOOST_CHECK_EQUAL(render(K, v[5]), "GTTCC");
    BOOST_CHECK_EQUAL(render(K, v[10]), "GGGTG");
    BOOST_CHECK_EQUAL(render(K, v[15]), "TCGCC");
    BOOST_CHECK_EQUAL(render(K, v[20]), "GCTGG");
    BOOST_CHECK_EQUAL(render(K, v[25]), "ATCGG");
    BOOST_CHECK_EQUAL(render(K, v[30]), "GAACC");
    BOOST_CHECK_EQUAL(render(K, v[35]), "TGGGC");
    BOOST_CHECK_EQUAL(render(K, v[40]), "GAGAC");
    BOOST_CHECK_EQUAL(render(K, v[45]), "AGTGG");
    BOOST_CHECK_EQUAL(render(K, v[50]), "AGCTG");
}

BOOST_AUTO_TEST_CASE( test_make_3 )
{
    using namespace std;
    using namespace varoom::kmers;

    const string s = "AGCCAGTTCCGGGTGTCGCCGCTGGATCGGNCCTGGAACCTGGGCGAGACAGTGGAGCTG";
    const size_t K = 5;
    vector<kmer> v, w;
    make(s, K, v, w);
    BOOST_CHECK_EQUAL(v.size(), 51);
    BOOST_CHECK_EQUAL(w.size(), 51);
    BOOST_CHECK_EQUAL(render(K, v[0]), "AGCCA");
    BOOST_CHECK_EQUAL(render(K, w[0]), "TGGCT");
    BOOST_CHECK_EQUAL(render(K, v[5]), "GTTCC");
    BOOST_CHECK_EQUAL(render(K, w[5]), "GGAAC");
    BOOST_CHECK_EQUAL(render(K, v[10]), "GGGTG");
    BOOST_CHECK_EQUAL(render(K, w[10]), "CACCC");
    BOOST_CHECK_EQUAL(render(K, v[15]), "TCGCC");
    BOOST_CHECK_EQUAL(render(K, w[15]), "GGCGA");
    BOOST_CHECK_EQUAL(render(K, v[20]), "GCTGG");
    BOOST_CHECK_EQUAL(render(K, w[20]), "CCAGC");
    BOOST_CHECK_EQUAL(render(K, v[25]), "ATCGG");
    BOOST_CHECK_EQUAL(render(K, w[25]), "CCGAT");
    BOOST_CHECK_EQUAL(render(K, v[30]), "GAACC");
    BOOST_CHECK_EQUAL(render(K, w[30]), "GGTTC");
    BOOST_CHECK_EQUAL(render(K, v[35]), "TGGGC");
    BOOST_CHECK_EQUAL(render(K, w[35]), "GCCCA");
    BOOST_CHECK_EQUAL(render(K, v[40]), "GAGAC");
    BOOST_CHECK_EQUAL(render(K, w[40]), "GTCTC");
    BOOST_CHECK_EQUAL(render(K, v[45]), "AGTGG");
    BOOST_CHECK_EQUAL(render(K, w[45]), "CCACT");
    BOOST_CHECK_EQUAL(render(K, v[50]), "AGCTG");
    BOOST_CHECK_EQUAL(render(K, w[50]), "CAGCT");
}

BOOST_AUTO_TEST_CASE( test_make_4 )
{
    using namespace std;
    using namespace varoom::kmers;

    const string s = "AGCCAGTTCCGGGTGTCGCCGCTGGATCGGNCCTGGAACCTGGGCGAGACAGTGGAGCTG";
    const size_t K = 5;
    vector<kmer_and_pos> v, w;
    make(s, K, v, w);
    BOOST_CHECK_EQUAL(v.size(), 51);
    BOOST_CHECK_EQUAL(w.size(), 51);
    BOOST_CHECK_EQUAL(render(K, v[0].first), "AGCCA");
    BOOST_CHECK_EQUAL(v[0].second, 0);
    BOOST_CHECK_EQUAL(render(K, w[0].first), "TGGCT");
    BOOST_CHECK_EQUAL(w[0].second, 55);
    BOOST_CHECK_EQUAL(render(K, v[5].first), "GTTCC");
    BOOST_CHECK_EQUAL(v[5].second, 5);
    BOOST_CHECK_EQUAL(render(K, w[5].first), "GGAAC");
    BOOST_CHECK_EQUAL(w[5].second, 50);
    BOOST_CHECK_EQUAL(render(K, v[45].first), "AGTGG");
    BOOST_CHECK_EQUAL(v[45].second, 50);
    BOOST_CHECK_EQUAL(render(K, w[45].first), "CCACT");
    BOOST_CHECK_EQUAL(w[45].second, 5);
    BOOST_CHECK_EQUAL(render(K, v[50].first), "AGCTG");
    BOOST_CHECK_EQUAL(v[50].second, 55);
    BOOST_CHECK_EQUAL(render(K, w[50].first), "CAGCT");
    BOOST_CHECK_EQUAL(w[50].second, 0);
}
