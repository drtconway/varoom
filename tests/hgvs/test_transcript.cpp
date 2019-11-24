#include "varoom/hgvs/transcript.hpp"

#include <unordered_map>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hgvs_positions tests
#include <boost/test/unit_test.hpp>

using namespace std;
using namespace varoom::hgvs;
using namespace nlohmann;

namespace // anonymous
{
    json kras = R"({
   "accession" : "NM_033360.3",
   "cdsEnd" : 25398318,
   "cdsStart" : 25368374,
   "chr" : "chr12",
   "exons" : [
      [25357722, 25362845],
      [25368370, 25368494],
      [25378547, 25378707],
      [25380167, 25380346],
      [25398207, 25398329],
      [25403684, 25403865]
   ],
   "strand" : "-",
   "txEnd" : 25403865,
   "txStart" : 25357722
})"_json;

    json msh2 = R"({
   "accession" : "NM_000251.2",
   "cdsEnd" : 47710088,
   "cdsStart" : 47630330,
   "chr" : "chr2",
   "exons" : [
      [47630205, 47630541],
      [47635539, 47635694],
      [47637232, 47637511],
      [47639552, 47639699],
      [47641407, 47641557],
      [47643434, 47643568],
      [47656880, 47657080],
      [47672686, 47672796],
      [47690169, 47690293],
      [47693796, 47693947],
      [47698103, 47698201],
      [47702163, 47702409],
      [47703505, 47703710],
      [47705410, 47705658],
      [47707834, 47708010],
      [47709917, 47710367]
   ],
   "strand" : "+",
   "txEnd" : 47710367,
   "txStart" : 47630205
}
)"_json;

    json brca2 = R"({
   "accession" : "NM_000059.3",
   "cdsEnd" : 32972907,
   "cdsStart" : 32890597,
   "chr" : "chr13",
   "exons" : [
      [32889616, 32889804],
      [32890558, 32890664],
      [32893213, 32893462],
      [32899212, 32899321],
      [32900237, 32900287],
      [32900378, 32900419],
      [32900635, 32900750],
      [32903579, 32903629],
      [32905055, 32905167],
      [32906408, 32907524],
      [32910401, 32915333],
      [32918694, 32918790],
      [32920963, 32921033],
      [32928997, 32929425],
      [32930564, 32930746],
      [32931878, 32932066],
      [32936659, 32936830],
      [32937315, 32937670],
      [32944538, 32944694],
      [32945092, 32945237],
      [32950806, 32950928],
      [32953453, 32953652],
      [32953886, 32954050],
      [32954143, 32954282],
      [32968825, 32969070],
      [32971034, 32971181],
      [32972298, 32973809]
   ],
   "strand" : "+",
   "txEnd" : 32973809,
   "txStart" : 32889616
})"_json;

    json pbx1 = R"({
   "accession" : "NM_001204961.1",
   "cdsEnd" : 164790820,
   "cdsStart" : 164529059,
   "chr" : "chr1",
   "exons" : [
      [164528596, 164529250],
      [164532474, 164532548],
      [164761730, 164761975],
      [164768935, 164769126],
      [164776778, 164776914],
      [164781226, 164781386],
      [164790773, 164790863],
      [164815820, 164821067]
   ],
   "strand" : "+",
   "txEnd" : 164821067,
   "txStart" : 164528596
})"_json;

}
// namespace anonymous

BOOST_AUTO_TEST_CASE( testTxPos1 )
{
    // cout << "MSH2" << endl;
    transcript tx = transcript::make(msh2);
    {
        // chr2:g.47630152C>T -> NM_000251.2:c.-179C>T   
        hgvsg_locus g(47630152);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR5);
        BOOST_CHECK_EQUAL(rp, -179);
        BOOST_CHECK_EQUAL(rr, 0);
    }
    {
        // chr2:g.47710431C>G -> NM_000251.2:c.*343C>G	
        hgvsg_locus g(47710431);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, 343);
        BOOST_CHECK_EQUAL(rr, 0);
    }
    {
        // chr2:g.47630221G>A -> NM_000251.2:c.-110G>A
        hgvsg_locus g(47630221);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR5);
        BOOST_CHECK_EQUAL(rp, -110);
        BOOST_CHECK_EQUAL(rr, 0);
    }
    {
        // chr2:g.47710151G>C -> NM_000251.2:c.*63G>C
        hgvsg_locus g(47710151);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, 63);
        BOOST_CHECK_EQUAL(rr, 0);
    }
    {
        hgvsg_locus g(47635539);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, INTRON);
        BOOST_CHECK_EQUAL(rp, 212);
        BOOST_CHECK_EQUAL(rr, -1);
    }
    {
        hgvsg_locus g(47641407);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, INTRON);
        BOOST_CHECK_EQUAL(rp, 793);
        BOOST_CHECK_EQUAL(rr, -1);
    }
    {
        hgvsg_locus g(47641558);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, INTRON);
        BOOST_CHECK_EQUAL(rp, 942);
        BOOST_CHECK_EQUAL(rr, 1);
    }
}
BOOST_AUTO_TEST_CASE( testTxPos2 )
{
    // cout << "BRCA2" << endl;
    transcript tx = transcript::make(brca2);
    {
        // chr13:g.32889968G>A -> NM_000059.3:c.-40+164G>A	
        hgvsg_locus g(32889968);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR5);
        BOOST_CHECK_EQUAL(rp, -40);
        BOOST_CHECK_EQUAL(rr, 164);
    }
    {
        // chr13:g.32890392C>T	NM_000059.3:c.-39-167C>T	
        hgvsg_locus g(32890392);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR5);
        BOOST_CHECK_EQUAL(rp, -39);
        BOOST_CHECK_EQUAL(rr, -167);
    }
}

BOOST_AUTO_TEST_CASE( testTxPos3 )
{
    // cout << "PBX1" << endl;
    transcript tx = transcript::make(pbx1);
    {
        // chr1:g.164790997G>T -> NM_001204961.1:c.*43+134G>T	
        hgvsg_locus g(164790997);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, 43);
        BOOST_CHECK_EQUAL(rr, 134);
    }
    {
        // chr1:g.164815744A>T	NM_001204961.1:c.*44-77A>T	
        hgvsg_locus g(164815744);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, 44);
        BOOST_CHECK_EQUAL(rr, -77);
    }
}

BOOST_AUTO_TEST_CASE( testTxNeg )
{
    // cout << "KRAS" << endl;
    transcript tx = transcript::make(kras);
    {
        // chr12:g.25403868del -> NM_033360.3:c.-195del
        hgvsg_locus g(25403868);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR5);
        BOOST_CHECK_EQUAL(rp, -195);
        BOOST_CHECK_EQUAL(rr, 0);
    }
    {
        // chr12:g.25357720del -> NM_033360.3:c.*5130del
        hgvsg_locus g(25357720);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, 5130);
        BOOST_CHECK_EQUAL(rr, 0);
    }
    {
        // chr12:g.25362777A>G -> NM_033360.3:c.*73T>C
        hgvsg_locus g(25362777);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, 73);
        BOOST_CHECK_EQUAL(rr, 0);
    }
    {
        // chr12:g.25368206C>T -> NM_033360.3:c.*4+165G>A	
        hgvsg_locus g(25368206);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, 4);
        BOOST_CHECK_EQUAL(rr, 165);
    }
    {
        // chr12:g.25363051A>C -> NM_033360.2:c.*5-206T>G	
        hgvsg_locus g(25363051);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, 5);
        BOOST_CHECK_EQUAL(rr, -206);
    }
    {
        // chr12:g.25398325C>T -> NM_033360.2:c.-7G>A
        hgvsg_locus g(25398325);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR5);
        BOOST_CHECK_EQUAL(rp, -7);
        BOOST_CHECK_EQUAL(rr, 0);
    }
    {
        // chr12:g.25398405C>T -> NM_033360.3:c.-11-76G>A	
        hgvsg_locus g(25398405);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR5);
        BOOST_CHECK_EQUAL(rp, -11);
        BOOST_CHECK_EQUAL(rr, -76);
    }
    {
        // chr12:g.25403680del -> NM_033360.3:c.-12+5del
        hgvsg_locus g(25403680);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR5);
        BOOST_CHECK_EQUAL(rp, -12);
        BOOST_CHECK_EQUAL(rr, 5);
    }
    {
        // chr12:g.25398309A>G -> NM_033360.3:c.10T>C
        hgvsg_locus g(25398309);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, CODING);
        BOOST_CHECK_EQUAL(rp, 10);
        BOOST_CHECK_EQUAL(rr, 0);
    }
    {
        // chr12:g.25398207C>G -> NM_033360.3:c.111+1G>C	
        hgvsg_locus g(25398207);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, INTRON);
        BOOST_CHECK_EQUAL(rp, 111);
        BOOST_CHECK_EQUAL(rr, 1);
    }
    {
        // chr12:g.25380360A>T -> NM_033360.3:c.112-14T>A	
        hgvsg_locus g(25380360);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, INTRON);
        BOOST_CHECK_EQUAL(rp, 112);
        BOOST_CHECK_EQUAL(rr, -14);
    }
}

BOOST_AUTO_TEST_CASE( testBoundariesPos )
{
    return;

    unordered_map<string,tx_zone> Z;
    Z["coding"] = CODING;
    Z["intron"] = INTRON;
    Z["utr3"] = UTR3;
    Z["utr5"] = UTR5;

    transcript msh2Tx = transcript::make(msh2);
    json msh2Tests = R"([
        [47635539, "intron", 212, -1],
        [47641407, "intron", 793, -1],
        [47641558, "intron", 942, 1],
        [47672686, "intron", 1277, -1]
    ])"_json;

    for (size_t i = 0; i < msh2Tests.size(); ++i)
    {
        uint64_t gj = msh2Tests[i][0];
        string zj = msh2Tests[i][1];
        int64_t pj = msh2Tests[i][2];
        int64_t rj = msh2Tests[i][3];

        hgvsg_locus g(gj);
        hgvsc_locus c = msh2Tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, Z[zj]);
        BOOST_CHECK_EQUAL(rp, pj);
        BOOST_CHECK_EQUAL(rr, rj);
    }

    transcript brca2Tx = transcript::make(brca2);
    json brca2Tests = R"([
        [32890665, "intron", 67, 1],
        [32893213, "intron", 68, -1],
        [32900288, "intron", 475, 1],
        [32972298, "intron", 9649, -1]
    ])"_json;

    for (size_t i = 0; i < brca2Tests.size(); ++i)
    {
        uint64_t gj = brca2Tests[i][0];
        string zj = brca2Tests[i][1];
        int64_t pj = brca2Tests[i][2];
        int64_t rj = brca2Tests[i][3];

        hgvsg_locus g(gj);
        hgvsc_locus c = brca2Tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, Z[zj]);
        BOOST_CHECK_EQUAL(rp, pj);
        BOOST_CHECK_EQUAL(rr, rj);
    }

    transcript pbx1Tx = transcript::make(pbx1);
    json pbx1Tests = R"([
        [164529251, "intron", 191, 1],
        [164532474, "intron", 192, -1],
        [164790864, "utr3", 43, 1],
        [164815820,  "utr3", 44, -1]
    ])"_json;

    for (size_t i = 0; i < pbx1Tests.size(); ++i)
    {
        uint64_t gj = pbx1Tests[i][0];
        string zj = pbx1Tests[i][1];
        int64_t pj = pbx1Tests[i][2];
        int64_t rj = pbx1Tests[i][3];

        hgvsg_locus g(gj);
        hgvsc_locus c = pbx1Tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, Z[zj]);
        BOOST_CHECK_EQUAL(rp, pj);
        BOOST_CHECK_EQUAL(rr, rj);
    }
}

