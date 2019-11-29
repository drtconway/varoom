#include "varoom/sam/pileup.hpp"

#include <iostream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE sam_pileup tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
    struct pile
    {
        std::string chr;
        std::uint32_t pos;
        std::string base;
        size_t count;
    };

    class pileup_tester : public varoom::sam_pileup
    {
    public:
        std::vector<pile> piles;

        virtual void output_pileup(const std::string& p_chr, const std::uint32_t& p_pos, const std::string& p_base, const size_t& p_count)
        {
            pile p;
            p.chr = p_chr;
            p.pos = p_pos;
            p.base = p_base;
            p.count = p_count;
            piles.push_back(p);
        }
    };
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( test1 )
{
    using namespace std;
    using namespace varoom;

    pileup_tester P;
    BOOST_CHECK_EQUAL(P.piles.size(), 0);

    string chr = "18";
    uint32_t pos = 16097;
    string seq = "TTTACACTTTGTGTATAGCAGGGAATCTGTGTCTAATTTGTAGTATTTCATGCTTCTAGGTTTTCATGGCAGTTGAGATGTAAGAATAACAATAATGTTG";
    string cigar = "100M";
    P.add_alignment(chr, pos, seq, cigar);
    BOOST_CHECK_EQUAL(P.piles.size(), 0);

    P.end();
    BOOST_CHECK_EQUAL(P.piles.size(), 100);
    BOOST_CHECK_EQUAL(P.piles[0].chr, "18");
    BOOST_CHECK_EQUAL(P.piles[0].pos, 16097);
    BOOST_CHECK_EQUAL(P.piles[0].base, "T");
    BOOST_CHECK_EQUAL(P.piles[0].count, 1);
    BOOST_CHECK_EQUAL(P.piles[5].chr, "18");
    BOOST_CHECK_EQUAL(P.piles[5].pos, 16102);
    BOOST_CHECK_EQUAL(P.piles[5].base, "A");
    BOOST_CHECK_EQUAL(P.piles[5].count, 1);
    BOOST_CHECK_EQUAL(P.piles[99].chr, "18");
    BOOST_CHECK_EQUAL(P.piles[99].pos, 16196);
    BOOST_CHECK_EQUAL(P.piles[99].base, "G");
    BOOST_CHECK_EQUAL(P.piles[99].count, 1);
}

BOOST_AUTO_TEST_CASE( test2Ins )
{
    using namespace std;
    using namespace varoom;

    pileup_tester P;
    BOOST_CHECK_EQUAL(P.piles.size(), 0);

    string chr = "18";
    uint32_t pos = 29769796;
    string seq = "AGCCAGGCACTTAATTGTTGACACAGGAAAAAAAAAAAAAAAGCCAAAAAACCAACTCCGCAGGGCGGTGTTAAA";
    string cigar = "27M1I47M";
    P.add_alignment(chr, pos, seq, cigar);
    BOOST_CHECK_EQUAL(P.piles.size(), 0);

    P.end();

    BOOST_CHECK_EQUAL(P.piles[26].pos, 29769822);
    BOOST_CHECK_EQUAL(P.piles[26].base, "G");
    BOOST_CHECK_EQUAL(P.piles[27].pos, 29769822);
    BOOST_CHECK_EQUAL(P.piles[27].base, "^A");
    BOOST_CHECK_EQUAL(P.piles[28].pos, 29769823);
    BOOST_CHECK_EQUAL(P.piles[28].base, "A");
    BOOST_CHECK_EQUAL(P.piles[29].pos, 29769824);
    BOOST_CHECK_EQUAL(P.piles[29].base, "A");
    BOOST_CHECK_EQUAL(P.piles[30].pos, 29769825);
    BOOST_CHECK_EQUAL(P.piles[30].base, "A");
}

BOOST_AUTO_TEST_CASE( test3Del )
{
    using namespace std;
    using namespace varoom;

    pileup_tester P;
    BOOST_CHECK_EQUAL(P.piles.size(), 0);

    string chr = "18";
    uint32_t pos = 29769954;
    string seq = "TAAATTATATGTTTGCCTTTAAAGAAGTCAATTCAACCCTGAATGTCATAGTTAGCTACTTTCAACTGGAAGCTA";
    string cigar = "8M1D67M";

    P.add_alignment(chr, pos, seq, cigar);
    BOOST_CHECK_EQUAL(P.piles.size(), 0);

    P.end();

    BOOST_CHECK_EQUAL(P.piles[7].pos, 29769961);
    BOOST_CHECK_EQUAL(P.piles[7].base, "*1");
    BOOST_CHECK_EQUAL(P.piles[8].pos, 29769961);
    BOOST_CHECK_EQUAL(P.piles[8].base, "T");
    BOOST_CHECK_EQUAL(P.piles[9].pos, 29769963);
    BOOST_CHECK_EQUAL(P.piles[9].base, "A");
}
