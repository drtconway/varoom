#include "varoom/vcf/vcf_info.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE vcf_info tests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( test1 )
{
    using namespace std;
    using namespace varoom::vcf;

    const std::string txt = "AF_ESP=0.00023;AF_EXAC=0.00018;ALLELEID=249324;CLNDISDB=MedGen:CN169374;CLNDN=not_specified;CLNHGVS=NC_000001.10:g.981328C>G;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=AGRN:375790;MC=SO:0001627|intron_variant;ORIGIN=1;RS=199785742";
    vcf_info ifo(txt);
    BOOST_CHECK_EQUAL(ifo.size(), 14);

    vector<string> ks;
    ifo.keys(ks);

    BOOST_CHECK_EQUAL(ks.size(), 14);
    BOOST_CHECK_EQUAL(ks[0], "AF_ESP");
    BOOST_CHECK_EQUAL(ks[1], "AF_EXAC");
    BOOST_CHECK_EQUAL(ks[2], "ALLELEID");
    BOOST_CHECK_EQUAL(ks[3], "CLNDISDB");
    BOOST_CHECK_EQUAL(ks[4], "CLNDN");
    BOOST_CHECK_EQUAL(ks[5], "CLNHGVS");
    BOOST_CHECK_EQUAL(ks[6], "CLNREVSTAT");
    BOOST_CHECK_EQUAL(ks[7], "CLNSIG");
    BOOST_CHECK_EQUAL(ks[8], "CLNVC");
    BOOST_CHECK_EQUAL(ks[9], "CLNVCSO");
    BOOST_CHECK_EQUAL(ks[10], "GENEINFO");
    BOOST_CHECK_EQUAL(ks[11], "MC");
    BOOST_CHECK_EQUAL(ks[12], "ORIGIN");
    BOOST_CHECK_EQUAL(ks[13], "RS");

    double af_esp;
    ifo.get("AF_ESP", af_esp);
    BOOST_CHECK_EQUAL(af_esp, 0.00023);

    uint64_t alleleid;
    ifo.get("ALLELEID", alleleid);
    BOOST_CHECK_EQUAL(alleleid, 249324);

    string geneinfo;
    ifo.get("CLNSIG", geneinfo);
    BOOST_CHECK_EQUAL(geneinfo, "Likely_benign");
}

