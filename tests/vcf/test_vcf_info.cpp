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

BOOST_AUTO_TEST_CASE( test2 )
{
    using namespace std;
    using namespace varoom::vcf;

    const std::string txt = "AC=1;AF=0.5;AN=2;BaseQRankSum=4.295;ClippingRankSum=-1.01;DB=true;FS=0.394;MLEAC=1;MLEAF=0.5;MQ=60;MQ0=0;MQRankSum=2.298;QD=15.02;ReadPosRankSum=-1.092;HGVSg=chr1:g.2488153A>G;HGVSc=NM_003820.2:c.50A>G;HGVSp=NP_003811.2:p.(Lys17Arg);gene=TNFRSF14";
    vcf_info ifo(txt);
    BOOST_CHECK_EQUAL(ifo.size(), 18);

    vector<string> ks;
    ifo.keys(ks);

    int64_t x;
    double y;
    string z;

    BOOST_CHECK_EQUAL(ks.size(), 18);

    BOOST_CHECK_EQUAL(ks[0],  "AC");
    BOOST_CHECK_EQUAL(ifo.get("AC", x), true);
    BOOST_CHECK_EQUAL(x, 1);

    BOOST_CHECK_EQUAL(ks[1],  "AF");
    BOOST_CHECK_EQUAL(ifo.get("AF", y), true);
    BOOST_CHECK_EQUAL(y, 0.5);

    BOOST_CHECK_EQUAL(ks[2],  "AN");
    BOOST_CHECK_EQUAL(ifo.get("AN", x), true);
    BOOST_CHECK_EQUAL(x, 2);

    BOOST_CHECK_EQUAL(ks[3],  "BaseQRankSum");
    BOOST_CHECK_EQUAL(ifo.get("BaseQRankSum", y), true);
    BOOST_CHECK_EQUAL(y, 4.295);

    BOOST_CHECK_EQUAL(ks[4],  "ClippingRankSum");
    BOOST_CHECK_EQUAL(ks[5],  "DB");
    BOOST_CHECK_EQUAL(ks[6],  "FS");
    BOOST_CHECK_EQUAL(ks[7],  "MLEAC");
    BOOST_CHECK_EQUAL(ks[8],  "MLEAF");
    BOOST_CHECK_EQUAL(ks[9],  "MQ");
    BOOST_CHECK_EQUAL(ks[10], "MQ0");
    BOOST_CHECK_EQUAL(ks[11], "MQRankSum");
    BOOST_CHECK_EQUAL(ks[12], "QD");
    BOOST_CHECK_EQUAL(ks[13], "ReadPosRankSum");
    BOOST_CHECK_EQUAL(ks[14], "HGVSg");
    BOOST_CHECK_EQUAL(ks[15], "HGVSc");

    BOOST_CHECK_EQUAL(ks[16], "HGVSp");
    BOOST_CHECK_EQUAL(ifo.get("HGVSp", z), true);
    BOOST_CHECK_EQUAL(z, "NP_003811.2:p.(Lys17Arg)");

    BOOST_CHECK_EQUAL(ks[17], "gene");
    BOOST_CHECK_EQUAL(ifo.get("gene", z), true);
    BOOST_CHECK_EQUAL(z, "TNFRSF14");
}

BOOST_AUTO_TEST_CASE( test3 )
{
    using namespace std;
    using namespace varoom::vcf;

    const std::string txt = "AC=1;AF=0.5;AN=2;BaseQRankSum=4.295;ClippingRankSum=-1.01;DB=true;FS=0.394;MLEAC=1;MLEAF=0.5;MQ=60;MQ0=0;MQRankSum=2.298;QD=15.02;ReadPosRankSum=-1.092;HGVSg=chr1:g.2488153A>G;HGVSc=NM_003820.2:c.50A>G;HGVSp=NP_003811.2:p.(Lys17Arg);gene=TNFRSF14";
    vcf_info ifo(txt);
    BOOST_CHECK_EQUAL(ifo.size(), 18);

    vector<string> ks;
    ifo.keys(ks);
    BOOST_CHECK_EQUAL(ks.size(), 18);

    string z;
    BOOST_CHECK_EQUAL(ifo.get("wibble", z), false);
}

BOOST_AUTO_TEST_CASE( test4 )
{
    using namespace std;
    using namespace varoom::vcf;

    const std::string fmt = "GT:AD:DP:GQ:PL:PMCAD:PMCADF:PMCADR:PMCBDIR:PMCDP:PMCFREQ:PMCRD:PMCRDF:PMCRDR";
    const std::string gtp = "0/1:159,177:336:99:5091,0,4233:196:103:93:Y:372:0.53:175:91:84";
    vcf_info ifo(fmt, gtp);
    BOOST_CHECK_EQUAL(ifo.size(), 14);

    string z;
    BOOST_CHECK_EQUAL(ifo.get("GT", z), true);
    BOOST_CHECK_EQUAL(z, "0/1");

    int64_t x;
    BOOST_CHECK_EQUAL(ifo.get("DP", x), true);
    BOOST_CHECK_EQUAL(x, 336);

    vector<int64_t> xs;
    BOOST_CHECK_EQUAL(ifo.get("AD", xs, ','), true);
    BOOST_CHECK_EQUAL(xs.size(), 2);
    BOOST_CHECK_EQUAL(xs[0], 159);
    BOOST_CHECK_EQUAL(xs[1], 177);
}
