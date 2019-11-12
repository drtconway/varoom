#include "varoom/vcf/vcf_info.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE vcf_info tests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( test1 )
{
    using namespace std;
    using namespace varoom::vcf;

    const std::string txt = "AF_ESP=0.00023;AF_EXAC=0.00018;ALLELEID=249324;CLNDISDB=MedGen:CN169374;CLNDN=not_specified;CLNHGVS=NC_000001.10:g.981328C>G;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=AGRN:375790;MC=SO:0001627|intron_variant;ORIGIN=1;RS=199785742";
    vcf_info_subtext ifo_mkr(txt);
    vcf_info ifo = ifo_mkr.make();
    BOOST_CHECK_EQUAL(ifo.size(), 14);

    size_t i = 0;
    BOOST_CHECK_EQUAL(ifo[i++].first, "AF_ESP");
    BOOST_CHECK_EQUAL(ifo[i++].first, "AF_EXAC");
    BOOST_CHECK_EQUAL(ifo[i++].first, "ALLELEID");
    BOOST_CHECK_EQUAL(ifo[i++].first, "CLNDISDB");
    BOOST_CHECK_EQUAL(ifo[i++].first, "CLNDN");
    BOOST_CHECK_EQUAL(ifo[i++].first, "CLNHGVS");
    BOOST_CHECK_EQUAL(ifo[i++].first, "CLNREVSTAT");
    BOOST_CHECK_EQUAL(ifo[i++].first, "CLNSIG");
    BOOST_CHECK_EQUAL(ifo[i++].first, "CLNVC");
    BOOST_CHECK_EQUAL(ifo[i++].first, "CLNVCSO");
    BOOST_CHECK_EQUAL(ifo[i++].first, "GENEINFO");
    BOOST_CHECK_EQUAL(ifo[i++].first, "MC");
    BOOST_CHECK_EQUAL(ifo[i++].first, "ORIGIN");
    BOOST_CHECK_EQUAL(ifo[i++].first, "RS");
    BOOST_CHECK_EQUAL(ifo.size(), i);
}

BOOST_AUTO_TEST_CASE( test4 )
{
    using namespace std;
    using namespace varoom::vcf;

    const std::string fmt = "GT:AD:DP:GQ:PL:PMCAD:PMCADF:PMCADR:PMCBDIR:PMCDP:PMCFREQ:PMCRD:PMCRDF:PMCRDR";
    const std::string gtp = "0/1:159,177:336:99:5091,0,4233:196:103:93:Y:372:0.53:175:91:84";

    vcf_info_subtext ifo_mkr(fmt, gtp);
    vcf_info ifo = ifo_mkr.make();
    BOOST_CHECK_EQUAL(ifo.size(), 14);
}
