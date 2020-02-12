#include "varoom/seq/genome_map.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE genome_map tests
#include <boost/test/unit_test.hpp>

using namespace std;
using namespace nlohmann;

namespace // anonymous
{
    json hg19 = R"({
        "chr1": 249250621,
        "chr2": 243199373,
        "chr3": 198022430,
        "chr4": 191154276,
        "chr5": 180915260,
        "chr6": 171115067,
        "chr7": 159138663,
        "chr8": 146364022,
        "chr9": 141213431,
        "chr10": 135534747,
        "chr11": 135006516,
        "chr12": 133851895,
        "chr13": 115169878,
        "chr14": 107349540,
        "chr15": 102531392,
        "chr16": 90354753,
        "chr17": 81195210,
        "chr18": 78077248,
        "chr19": 59128983,
        "chr20": 63025520,
        "chr21": 48129895,
        "chr22": 51304566,
        "chrX": 155270560,
        "chrY": 59373566,
        "chrM": 16571,

        "chr4_ctg9_hap1": 590426,
        "chr6_apd_hap1": 4622290,
        "chr6_cox_hap2": 4795371,
        "chr6_dbb_hap3": 4610396,
        "chr6_mann_hap4": 4683263,
        "chr6_mcf_hap5": 4833398,
        "chr6_qbl_hap6": 4611984,
        "chr6_ssto_hap7": 4928567,
        "chr17_ctg5_hap1": 1680828,

        "chr1_gl000191_random": 106433,
        "chr1_gl000192_random": 547496,
        "chr4_gl000193_random": 189789,
        "chr4_gl000194_random": 191469,
        "chr7_gl000195_random": 182896,
        "chr8_gl000196_random": 38914,
        "chr8_gl000197_random": 37175,
        "chr9_gl000198_random": 90085,
        "chr9_gl000199_random": 169874,
        "chr9_gl000200_random": 187035,
        "chr9_gl000201_random": 36148,
        "chr11_gl000202_random": 40103,
        "chr17_gl000203_random": 37498,
        "chr17_gl000204_random": 81310,
        "chr17_gl000205_random": 174588,
        "chr17_gl000206_random": 41001,
        "chr18_gl000207_random": 4262,
        "chr19_gl000208_random": 92689,
        "chr19_gl000209_random": 159169,
        "chr21_gl000210_random": 27682,
        "chrUn_gl000211": 166566,
        "chrUn_gl000212": 186858,
        "chrUn_gl000213": 164239,
        "chrUn_gl000216": 172294,
        "chrUn_gl000217": 172149,
        "chrUn_gl000218": 161147,
        "chrUn_gl000219": 179198,
        "chrUn_gl000220": 161802,
        "chrUn_gl000221": 155397,
        "chrUn_gl000222": 186861,
        "chrUn_gl000223": 180455,
        "chrUn_gl000224": 179693,
        "chrUn_gl000225": 211173,
        "chrUn_gl000226": 15008,
        "chrUn_gl000227": 128374,
        "chrUn_gl000228": 129120,
        "chrUn_gl000229": 19913,
        "chrUn_gl000230": 43691,
        "chrUn_gl000231": 27386,
        "chrUn_gl000232": 40652,
        "chrUn_gl000233": 45941,
        "chrUn_gl000234": 40531,
        "chrUn_gl000235": 34474,
        "chrUn_gl000236": 41934,
        "chrUn_gl000237": 45867,
        "chrUn_gl000238": 39939,
        "chrUn_gl000239": 33824,
        "chrUn_gl000240": 41933,
        "chrUn_gl000241": 42152,
        "chrUn_gl000242": 43523,
        "chrUn_gl000243": 43341,
        "chrUn_gl000244": 39929,
        "chrUn_gl000245": 36651,
        "chrUn_gl000246": 38154,
        "chrUn_gl000247": 36422,
        "chrUn_gl000248": 39786,
        "chrUn_gl000249": 38502
    })"_json;
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( hg19_1 )
{
    using namespace varoom::seq;

    genome_map m(hg19);

    {
        genome_map::chrom ch("chr1");
        size_t pos = 0;
        size_t p = m.chrom2genome(ch, pos);
        std::pair<genome_map::chrom,size_t> r = m.genome2chrom(p);
        BOOST_CHECK_EQUAL(r.first, ch);
        BOOST_CHECK_EQUAL(r.second, pos);
    }
    {
        genome_map::chrom ch("chr1");
        size_t pos = 123456;
        size_t p = m.chrom2genome(ch, pos);
        std::pair<genome_map::chrom,size_t> r = m.genome2chrom(p);
        BOOST_CHECK_EQUAL(r.first, ch);
        BOOST_CHECK_EQUAL(r.second, pos);
    }
    {
        genome_map::chrom ch("chr7");
        size_t pos = 0;
        size_t p = m.chrom2genome(ch, pos);
        std::pair<genome_map::chrom,size_t> r = m.genome2chrom(p);
        BOOST_CHECK_EQUAL(r.first, ch);
        BOOST_CHECK_EQUAL(r.second, pos);
    }
    {
        genome_map::chrom ch("chr7");
        size_t pos = 123456;
        size_t p = m.chrom2genome(ch, pos);
        std::pair<genome_map::chrom,size_t> r = m.genome2chrom(p);
        BOOST_CHECK_EQUAL(r.first, ch);
        BOOST_CHECK_EQUAL(r.second, pos);
    }
}

BOOST_AUTO_TEST_CASE( hg19_2 )
{
    const varoom::seq::genome_map& g = varoom::seq::genome_map::hg19();
    BOOST_CHECK_EQUAL(g.size(), 3136851001ULL);
}

BOOST_AUTO_TEST_CASE( hg38_1 )
{
    const varoom::seq::genome_map& g = varoom::seq::genome_map::hg38();
    BOOST_CHECK_EQUAL(g.size(), 3209286105ULL);
}