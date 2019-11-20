#include "varoom/hgvs/transcript.hpp"

#include <nlohmann/json.hpp>

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

    transcript make(const json& x)
    {
        tx_strand s = POS;
        if (x["strand"] == "-")
        {
            s = NEG;
        }
        string acc = x["accession"];
        string chr = x["chr"];
        uint64_t tss = x["txStart"];
        uint64_t tse = x["txEnd"];
        uint64_t cdss = x["cdsStart"];
        uint64_t cdse = x["cdsEnd"];

        genomic_locus txStart(tss);
        genomic_locus txEnd(tse);
        genomic_locus cdsStart(cdss);
        genomic_locus cdsEnd(cdse);

        vector<genomic_exon> exons;
        for (size_t i = 0; i < x["exons"].size(); ++i)
        {
            const json& y = x["exons"][i];
            uint64_t exs = y[0];
            uint64_t exe = y[1];
            exons.push_back(genomic_exon(genomic_locus(exs), genomic_locus(exe)));
        }
        return transcript(acc, chr, s, txStart, txEnd, cdsStart, cdsEnd, exons);
    }

}
// namespace anonymous

BOOST_AUTO_TEST_CASE( testTxPos1 )
{
    // cout << "MSH2" << endl;
    transcript tx = make(msh2);
    {
        // chr2:g.47630152C>T -> NM_000251.2:c.-179C>T   
        genomic_locus g0(47630152);
        hgvsg_locus g(static_cast<uint64_t>(g0) + 1);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR5);
        BOOST_CHECK_EQUAL(rp, -179);
        BOOST_CHECK_EQUAL(rr, 0);
    }
    {
        // chr2:g.47710431C>G -> NM_000251.2:c.*343C>G	
        genomic_locus g0(47710431);
        hgvsg_locus g(static_cast<uint64_t>(g0) + 1);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, 343);
        BOOST_CHECK_EQUAL(rr, 0);
    }
    {
        // chr2:g.47630221G>A -> NM_000251.2:c.-110G>A
        genomic_locus g0(47630221);
        hgvsg_locus g(static_cast<uint64_t>(g0) + 1);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR5);
        BOOST_CHECK_EQUAL(rp, -110);
        BOOST_CHECK_EQUAL(rr, 0);
    }
}
BOOST_AUTO_TEST_CASE( testTxPos2 )
{
    // cout << "BRCA2" << endl;
    transcript tx = make(brca2);
    {
        // chr13:g.32889968G>A -> NM_000059.3:c.-40+164G>A	
        genomic_locus g0(32889968);
        hgvsg_locus g(static_cast<uint64_t>(g0) + 1);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR5);
        BOOST_CHECK_EQUAL(rp, -40);
        BOOST_CHECK_EQUAL(rr, 164);
    }
    {
        // chr13:g.32890392C>T	NM_000059.3:c.-39-167C>T	
        genomic_locus g0(32890392);
        hgvsg_locus g(static_cast<uint64_t>(g0) + 1);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR5);
        BOOST_CHECK_EQUAL(rp, -39);
        BOOST_CHECK_EQUAL(rr, -167);
    }
}

BOOST_AUTO_TEST_CASE( testTxNeg )
{
    // cout << "KRAS" << endl;
    transcript tx = make(kras);
    {
        // chr12:g.25403868del -> NM_033360.3:c.-195del
        genomic_locus g0(25403868);
        hgvsg_locus g(static_cast<uint64_t>(g0) + 1);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR5);
        BOOST_CHECK_EQUAL(rp, -195);
        BOOST_CHECK_EQUAL(rr, 0);
    }
    {
        // chr12:g.25357720del -> NM_033360.3:c.*5130del
        genomic_locus g0(25357720);
        hgvsg_locus g(static_cast<uint64_t>(g0) + 1);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, 5130);
        BOOST_CHECK_EQUAL(rr, 0);
    }
    {
        // chr12:g.25362777A>G -> NM_033360.3:c.*73T>C
        genomic_locus g0(25362777);
        hgvsg_locus g(static_cast<uint64_t>(g0) + 1);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, 73);
        BOOST_CHECK_EQUAL(rr, 0);
    }
    {
        // chr12:g.25368206C>T -> NM_033360.3:c.*4+165G>A	
        genomic_locus g0(25368206);
        hgvsg_locus g(static_cast<uint64_t>(g0) + 1);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, 4);
        BOOST_CHECK_EQUAL(rr, 165);
    }
    {
        // chr12:g.25363051A>C -> NM_033360.2:c.*5-206T>G	
        genomic_locus g0(25363051);
        hgvsg_locus g(static_cast<uint64_t>(g0) + 1);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, 5);
        BOOST_CHECK_EQUAL(rr, -206);
    }
    if (0) {
        // chr12:g.25398405C>T -> NM_033360.3:c.-11-76G>A	
        genomic_locus g0(25398405);
        hgvsg_locus g(static_cast<uint64_t>(g0) + 1);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, -11);
        BOOST_CHECK_EQUAL(rr, -76);
    }
    if (0) {
        // chr12:g.25403680del -> NM_033360.3:c.-12+5del
        genomic_locus g0(25403680);
        hgvsg_locus g(static_cast<uint64_t>(g0) + 1);
        hgvsc_locus c = tx.to_hgvsc_locus(g);
        int64_t rp = static_cast<int64_t>(c.pos);
        int64_t rr = static_cast<int64_t>(c.rel);
        BOOST_CHECK_EQUAL(c.zone, UTR3);
        BOOST_CHECK_EQUAL(rp, -12);
        BOOST_CHECK_EQUAL(rr, 5);
    }
}
