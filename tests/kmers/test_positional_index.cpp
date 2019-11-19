#include "varoom/kmers/positional_index.hpp"

#include <iostream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE positional_index tests
#include <boost/test/unit_test.hpp>

using namespace std;
using namespace varoom;

namespace // anonymous
{

    const string seq =
        "AGTCTGTTTTTTTATTTCTACTTTTCCACAGTCTATTAGAAGCCTCTTCCCTTTATTGTCAAGGTTTTTTATTATTTCCA"
        "CAGGCTGTGGAAAACTTTAGTTAACACTGTGAATTACCTTTCCACAACTTGTGGGTAACTATAACTATTCTTTCAGGTTT"
        "TGTGGAAAACTGACTGGATTTGTGTTAAAATAGTCCTAGAATTATCCACAAGAAGGAACCTAGTATGACTGAAAATGAAC"
        "AAATTTTTTGGAATCGGGTCTTGGAATTAGCTCAGAGCCAATTAAAACAGGCAACGTATGAATTTTTTGTGCATGATGCA"
        "CGCTTATTAAAAGTTGAGAATCATGTGGCAACGATTTACTTAGATCAAATGAAAGAACTCTTTTGGGAAAAAAACCTTAA"
        "AGATGTTATTCTGACAGCTGGTTTTGAAGTTTATAATGCTCAAATTGCCGTTGACTATGTTTTTGAAGAAGATCTGATAA"
        "TTGACCAAAATCAAATCCCAAATAGCCAAAGCTACAATCAGCAAGTAATAACTCCTTTACCTGCTGTTACCTCAGACCTA"
        "AATCCCAAATATAGTTTTGAAAACTTTATTCAGGGTGATGAAAATCGTTGGGCTGTAGCTGCTTCTATAGCAGTAGCTAA"
        "TACACCAGGAACAACCTATAACCCTTTGTTTATCTGGGGAGGGCCTGGACTAGGAAAAACCCATTTGTTAAATGCCATTG"
        "GAAACTCTGTGCTATTAGAAAATCCCAATGCCCGTATAAAGTACATCACCGCTGAAAATTTCATTAATGAATTTGTGATT"
        "CATATTCGTTTAGACACTATGGACGAATTAAAAGAAAAGTTCCGTAATCTCGATTTACTGCTTATTGATGATATCCAATC"
        "GCTGGCGAAGAAAACATTATCTGGAACACAAGAAGAGTTCTTTAATACTTTTAATGCTCTTCATAACAATAATAAACAAA"
        "TCGTTCTAACCAGTGACCGCACACCAGATCACCTTAATGATTTAGAAGATCGATTGGTAACCCGTTTCAAATGGGGCTTA"
        "ACCGTTAATATTACGCCTCCAGATTTTGAAACTCGAGTAGCTATCTTAACGAATAAAATTCAGGAATATAATTTTATTTT"
        "CCCTCAAGATACTATTGAATATCTAGCTGGCCAATTCGATTCCAATGTGAGAGATTTGGAAGGCGCCTTAAAAGACATTA"
        "GTTTGGTTGCTAATTTTAAGCAAATTGATACCATTACGGTTGACGTTGCTGCAGAAGCTATTCGTGCCAGAAAACAAGAT"
        "GGACCTAAAATGACGGTCATTCCAATTGAAGAAATTCAAACGCAGGTTGGAAAATTCTATGACGTCACTGTCAAAGAAAT"
        "TAAGGCAACTAAACGTACGCAAGATATTGTGTTAGCAAGACAAGTCGCTATGTTTTTAGCACGTGAAATGACAGATAATA"
        "GCCTTCCCAAAATTGGTAAGGAATTTGGTGGTAGGGATCATTCCACTGTACTTCATGCCTACAATAAAATCAAAAACATG"
        "ATTAGCCAGGATGAAAGCCTTCGTATTGAAATTGAAACCATCAAAAATAAAATTAAGTAGCTTGTGGACAAGTTCTATTT"
        "TTAGTGACGAGTTATCCACAAGTTGTGAACAGTCTTCTTTCCTTATCCCTACTAGATAAATCAGACTTATCCACGTCATA"
        "CACAAGACCTACTACTACTACTAATTATTATACTTATCAATAAAGGAGTCCTCATGATTCAATTTTCCATAAATCGTACC"
        "CTTTTCATTCAAGCTTTAAATGCCACTAAACGTGGTATTAGCAGTAAAAATGCCATTCCTGTTCTTTCTACCATTAAGAT"
        "TAACGTTAGTTCATCTGATATCACTTTAACTGGTTCAAATGGACAAATTTCAATTGAAAATACCATTCCTGTAAGCAATG"
        "AAAATGCTGGACTATTAATCACATCTCCAGGGTCTATTCTTCTGGAAGCAAATTTCTTTATCAATATTATTTCTAGTTTG"
        "CCAGATGTTAGTTTGGATTTTAAAGAAATTGAACAACATCAAGTTGTTTTAACCAGTGGTAAATCAGAAATTACCTTAAA"
        "AGGAAAAGATGTTGATCAATACCCTCGGTTACAAGAAGTATCAACAGAAAATCCTTTGATCTTAAAAACAAAATTATTAA"
        "AGTCTATTATTGCTGAAACAGCTTTTGCAGCCAGTTTACAAGAAAGTCGTCCTATCTTAACAGGAGTTCATATTGTATTA"
        "AGTAACCATAAAGATTTTAAAGCCGTAGCAACTGACTCTCATCGTATGAGTCAACGTTTAATCACTTTGGATAATACTTC"
        "AGCAGATTTTGATGTGGTTATTCCAAGTAAATCTTTGAGAGAATTTTCAGCAGTATTTACAGATGATATTGAGACTGTTG"
        "AGGTATTTTTCTCACCAAGCCAAATCTTGTTCAGAAGTGAACATATTTCTTTCTATACACGTCTCTTAGAAGGAAATTAT"
        "CCCGATACAGACCGTTTGTTAATGACACAATTTGAGACAGAGGTTGTTTTCAATACCCAATCTCTTCGCCACGCTATGGA"
        "ACGTGCTTTCTTGATTTCGAACGCTACTCAAAATGGTACCGTTAAACTTGAAATTGCTCAAAATCATATTTCAGCTCATG"
        "TTAACTCACCGGAAGTTGGTAAAGTAAACGAAGATTTGGATATTGTTAGTCAATCTGGTAGTGATTTAACTATTAGTTTC"
        "AACCCAACTTACCTTATTGAATCTCTCAAGGCTATCAAGAGTGAAACAGTTAAAATTCATTTCTTGTCACCAGTACGACC"
        "ATTTACCTTGACACCAGGTGATGATGAAGAAAGCTTTATCCAGTTAATCACACCAGTCCGTACCAACTAAAAAGAAAAGG"
        "CTCCCTTTTAGGAGCTTTTTTTGTTATCATAAATGATGAAGATAATAAGAGTGAGGAAAAAAGATGTATCAAATTGGATC"
        "ACTTGTTGAAATGAAAAAACCGCACGCCTGTGTGATTAAAGAGACTGGTAAAAAATCTAATCAATGGAAAGTGCTTAGAG"
        "TAGGAGCTGATATTAAAATTCAATGCACTAACTGTCAGCACATCATTATGATGAGCCGTTACGACTTTGACCGAAAACTA"
        "AAAAAAGTCCTGCAACTTTAGAAATATTGATTTAGTAGGCTTTCTTATACATCTTGCAACCAATACTTGCCTAAATAATT"
        "GTTAGTATGCCTTTGGAAAATCAGGTATTCTAATGTTATCGAAAGAAGAAAGGTGGTCATAGAAAATGACAAAAGTTGCA"
        "GAACAATTAAAGCAATTACGAGTGAAACATCAATTATCTCAAGATGCTCTGGCAGAACAGTTATTTATTTCTCGGCAAGC"
        "CATATCAAAATGAGAAAATGGAGATACAATACCAGATTTGGATAATTTGGTCAGGTTAACTGAAATTTTTGACGTGAGCT"
        "TAGATGAGCTTGTTTTAGCTAAACCACATGAAGTTAAAGTTGAACGCATTTATGAAAACAAACCGCTTGATCTACAAAAA"
        "TACAATAAGCTCTATTGGTTTATTTTTCGAAATATTATTCTGTCTCTACTAATTATTTTAGCTATATTAACTATCTTAGA"
        "AGTTTTAGGGATACCTTTTGTTTCTAATTGGTTAATTTAAAGAAAAGTTGAAGAGTAATATGATGCTAAGCGAACCTGAA"
        "ATTTTCTTATCATCATACCGCTTTTTTCGTTTATTTTCTGTTATAATAGTTGTGATTGAAATTTTGAATGGAGACTTATT"
        "AAAATGGCTTTAACAGCAGGTATTGTGGGCTTACCTAATGTTGGTAAATCAACTTTATTTAATGCAATTACAAAAGCAGG"
        "GGCAGAAGCTGCTAATTATCCTTTTGCAACGATTGATCCTAATGTTGGGATGGTAGAGGTACCAGATGAACGTCTGCAAA"
        "AATTGACAGAGTTGATTACGCCTAAAAAAACCGTTCCAACAACCTTTGAGTTTACTGATATTGCGGGTATTGTTAAAGGA"
        "GCTTCTAAAGGAGAAGGGTTAGGTAATAAATTCTTGGCCAATATCCGTGAAGTAGATGCTATCGTACATGTCGTTCGTGC"
        "TTTTGATGATGAAAATGTTATGCGTGAACAAGGTCGTGAGGATGCTTTCGTGGATCCAATGGCTGATATTGATACTATCA"
        "ACCTTGAATTGATTTTGGCTGATTTAGAGTCTATTAATAAGCGTTATGCGCGTGTGGAAAAAATGGCTCGTACCCAAAAA"
        "GATAAGGATTCCGTGGCAGAATTTGCCGTTCTTGAAAAAATCAAACCTGTCTTAGAAGATGGTAAATCTGCTCGGACAAT";

}
// namespace anonymous

BOOST_AUTO_TEST_CASE( cons1 )
{
    typedef positional_index::seq_num_and_pos seq_num_and_pos;

    const size_t K = 25;
    positional_index idx(K);
    idx.add(seq);
    const string r = "CCCTCAAGATACTATTGAATATCTAGCTGGCCAATTCGATTCCAATGTGAGAGATTTGGAAGGCGCCTTAAAAGACATTA";
    vector<kmer_and_pos> xs;
    vector<kmer_and_pos> ys;
    kmers::make(r, K, xs, ys);
    for (size_t i = 0; i < xs.size(); ++i)
    {
        kmer x = xs[i].first;
        int32_t p = xs[i].second;
        vector<seq_num_and_pos> h;
        //cerr << kmers::render(K, x) << '\t' << p << endl;
        bool v = idx.hits(x, h);
        BOOST_CHECK_EQUAL(v, true);
        BOOST_REQUIRE_EQUAL(h.size(), 1);
        BOOST_CHECK_EQUAL(h[0].first, 0);
        BOOST_CHECK_EQUAL(h[0].second, 1120 + p);
    }
}

BOOST_AUTO_TEST_CASE( cons2 )
{
    typedef positional_index::seq_num seq_num;
    typedef positional_index::seq_pos seq_pos;

    const size_t K = 25;
    positional_index idx(K);
    idx.add(seq);
    const string r = "CCCTCAAGATACTATTGAATATCTAGCTGGCCAATTCGATTCCAATGTGAGAGATTTGGAAGGCGCCTTAAAAGACATTA";
    vector<kmer_and_pos> xs;
    vector<kmer_and_pos> ys;
    kmers::make(r, K, xs, ys);

    std::unordered_map<seq_num, std::unordered_map<seq_pos,size_t>> hx;
    bool v = idx.hits(xs, hx);
    BOOST_CHECK_EQUAL(v, true);
    BOOST_CHECK_EQUAL(hx.size(), 1);
    BOOST_CHECK_EQUAL(hx.count(0), 1);
    BOOST_REQUIRE_EQUAL(hx[0].size(), 1);
    BOOST_CHECK_EQUAL(hx[0].begin()->first, 1120);
    BOOST_CHECK_EQUAL(hx[0].begin()->second, 56);
}
