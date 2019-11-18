#include "varoom/seq/fastq_pair.hpp"

#include "varoom/util/files.hpp"

#include <initializer_list>
#include <fstream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE fastq_pair tests
#include <boost/test/unit_test.hpp>

using namespace std;
using namespace varoom;

namespace // anonymous
{
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( test_fastq_pair )
{
    using namespace std;
    using namespace varoom::seq;

    input_file_holder_ptr inp1 = files::in("tests/data/SRR2912526-100_1.fastq.gz");
    fastq_reader r1(**inp1);
    input_file_holder_ptr inp2 = files::in("tests/data/SRR2912526-100_2.fastq.gz");
    fastq_reader r2(**inp2);
    size_t n = 0;
    for (fastq_pair<fastq_reader,fastq_reader> r1r2(r1, r2); r1r2.more(); ++r1r2, ++n)
    {
        pair<const fastq_read&,const fastq_read&> rp = *r1r2;
        switch (n)
        {
            case 0:
            {
                const string s1 = "GGATAAAGTTTGGTAACATTGTGGATTATTTTTCACAGCTTGTGGAAAATTCTTGCTATCTATGGTAAAATATCTCTAGTATTAAACTTTTAAATAGAAAAGGAGGAGAAAGGATTGAAAGAAAAACAATTTTGGAATCGTATATTAGAATTAGCACAAGAAAGACTGACTCGATCCATGTATGATTTCTATGCTATTCAAGCTGAACTCATCAAGGTAGAGGAAAATGTTGCCACTATATTTCTACCTCGCTCTGAAATGGAAATGGTNTGGNANNAACAANNANNNNNNNNNNNNNNNN";
                const string s2 = "NNNNCATNNNNNNNNNTTNNNNNNNNNNNNNNNNNNNCANNNAANNNNCGNCCAGNNNCNGNNNNNNNANGTNNNCNNNNNNNNNANNNTNNNNNNNNNNNNNGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
                BOOST_CHECK_EQUAL(std::get<1>(rp.first), s1);
                BOOST_CHECK_EQUAL(std::get<1>(rp.second), s2);
                break;
            }
            case 24:
            {
                const string s1 = "CGATAGGTCTTTTTAAACTTTTCCATTTCCCCAAGTCTTAGGTGATCAAGAAAGTCATTAATAAAGCTTTCGGCAGGGATATATTTAACACGCGCATTAGGAATATTTTTTAGAATTTCATTTCCAATAGCGTTTAATAAGTGAGTCTTACCAAGGCCTGGTCCTCCATAGATAAAAAGAGGGTTATAGGTCAGAGCCAAATCTTCAGAGACAGCTAAAGCGGCTGATACAGCCCAAACATTTCCATCCCCTTGAATAAAGTTATCAAAGGTATACTTTTCTTTTAATCCCGTATCTGAAT";
                const string s2 = "GTAAAGGAGGAGAAAGGATTGAAAGAAAAACAATTTTGGAATCGTATATTAGAATTTGCACAAGAAAGACTGACTCGATCCATGTATGATTTCTATGCTATTCAAGCTGAACTCATCAAGGTAGAGGAAAATGTTGCCACTATATTTCTACCTCGCTCTGAAATGGAAATGGTCTGGGAAAAACAACTAAAAGATATTATT";
                BOOST_CHECK_EQUAL(std::get<1>(rp.first), s1);
                BOOST_CHECK_EQUAL(std::get<1>(rp.second), s2);
                break;
            }
            default:
            {
                // do nothing
            }
        }
    }
    BOOST_CHECK_EQUAL(n, 25);
}

