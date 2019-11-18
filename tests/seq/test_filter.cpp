#include "varoom/seq/fastq_filter.hpp"

#include "varoom/kmers.hpp"
#include "varoom/seq/fastq.hpp"
#include "varoom/seq/fastq_pair.hpp"
#include "varoom/util/files.hpp"

#include <iostream>
#include <unordered_set>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE filter tests
#include <boost/test/unit_test.hpp>

using namespace std;
using namespace varoom;
using namespace varoom::seq;

namespace // anonymous
{
    class catcher
    {
    public:
        catcher(const string& p_seq, size_t p_k)
            : m_k(p_k)
        {
            vector<kmer> xs;
            kmers::make(p_seq, p_k, xs);
            m_kmers.insert(xs.begin(), xs.end());
        }

        bool overlaps(const fastq_read& p_r) const
        {
            vector<kmer> xs;
            kmers::make(std::get<1>(p_r), m_k, xs);
            for (auto itr = xs.begin(); itr != xs.end(); ++itr)
            {
                if (m_kmers.count(*itr))
                {
                    return true;
                }
            }
            return false;
        }

    private:
        const size_t m_k;
        unordered_set<kmer> m_kmers;
    };

    vector<string> res{
        "SRR2912526.1 1 length=301",
        "SRR2912526.2 2 length=301",
        "SRR2912526.3 3 length=301",
        "SRR2912526.4 4 length=301",
        "SRR2912526.5 5 length=301",
        "SRR2912526.8 8 length=301",
        "SRR2912526.10 10 length=301",
        "SRR2912526.11 11 length=301",
        "SRR2912526.13 13 length=301",
        "SRR2912526.14 14 length=301",
        "SRR2912526.15 15 length=301",
        "SRR2912526.16 16 length=301",
        "SRR2912526.17 17 length=301",
        "SRR2912526.19 19 length=301",
        "SRR2912526.20 20 length=301",
        "SRR2912526.25 25 length=301"
    };
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( test_single )
{
    using namespace std;
    using namespace varoom::seq;

    const size_t K = 8;
    catcher cx("TCCCCTTGAATAAAGTTATCAAAGGTATACTTTTCTTTT", K);

    input_file_holder_ptr inp1 = files::in("tests/data/SRR2912526-100_1.fastq.gz");
    fastq_reader r1(**inp1);
    fastq_filter<fastq_reader> s1(r1, [cx](const fastq_read& rd) { return cx.overlaps(rd); });
    size_t n = 0;
    while (s1.more())
    {
        BOOST_CHECK_EQUAL(std::get<0>(*s1), res[n]);
        ++s1;
        ++n;
    }
    BOOST_CHECK_EQUAL(n, res.size());
}

BOOST_AUTO_TEST_CASE( test_pair )
{
    using namespace std;
    using namespace varoom::seq;

    typedef fastq_pair<fastq_reader,fastq_reader> fastq_pair_reader;
    typedef fastq_pair_reader::item_type read_pair;

    const size_t K = 8;
    catcher cx("TCCCCTTGAATAAAGTTATCAAAGGTATACTTTTCTTTT", K);

    input_file_holder_ptr inp1 = files::in("tests/data/SRR2912526-100_1.fastq.gz");
    fastq_reader r1(**inp1);
    input_file_holder_ptr inp2 = files::in("tests/data/SRR2912526-100_2.fastq.gz");
    fastq_reader r2(**inp2);
    fastq_pair_reader r1r2(r1, r2);

    std::function<bool(const read_pair&)> f = [cx](const read_pair& rd) { return cx.overlaps(rd.first); };

    fastq_filter<fastq_pair_reader> s1(r1r2, f);

    size_t n = 0;
    while (s1.more())
    {
        BOOST_CHECK_EQUAL(std::get<0>((*s1).first), res[n]);
        ++s1;
        ++n;
    }
    BOOST_CHECK_EQUAL(n, res.size());
}
