#include "varoom/seq/locus_stream.hpp"

#include <iostream>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE locus_stream tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
    class test_locus_stream : public varoom::seq::locus_stream<int>
    {
    public:
        test_locus_stream()
            : locus_stream<int>(0), m_chrs{"1", "2", "10", "20"}, m_cur_chr(0), m_cur_pos(0)
        {
            next();
        }

    private:
        virtual void next()
        {
            ++m_cur_pos;
            ++m_data;
            
            if (m_cur_pos < 5)
            {
                m_locus = varoom::seq::locus_id(m_chrs[m_cur_chr], m_cur_pos);
                return;
            }
            ++m_cur_chr;
            m_cur_pos = 1;
            if (m_cur_chr == m_chrs.size())
            {
                m_more = false;
                return;
            }
            m_locus = varoom::seq::locus_id(m_chrs[m_cur_chr], m_cur_pos);
        }

        const std::vector<std::string> m_chrs;
        size_t m_cur_chr;
        size_t m_cur_pos;
    };
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( test_locus_id )
{
    using namespace varoom::seq;

    locus_id a("1", 1234);
    locus_id b("2", 2345);
    locus_id c("2", 3456);
    locus_id d("10", 4567);
    locus_id e("234", 5678);
    locus_id f("1000", 6789);

    BOOST_CHECK_EQUAL(a < b, true);
    BOOST_CHECK_EQUAL(b < c, true);
    BOOST_CHECK_EQUAL(c < d, true);
    BOOST_CHECK_EQUAL(d < e, true);
    BOOST_CHECK_EQUAL(e < f, false);

    std::vector<locus_id> xs{f, e, d, c, b, a};
    std::sort(xs.begin(), xs.end());

    BOOST_CHECK_EQUAL(xs[0] == a, true);
}

BOOST_AUTO_TEST_CASE( test_locus_stream_1 )
{
    using namespace std;
    bool first = true;
    varoom::seq::locus_id p;
    int i = 1;
    for (test_locus_stream r; r.more(); ++r, ++i)
    {
        BOOST_CHECK_EQUAL(r.data(), i);
        if (!first)
        {
            BOOST_CHECK_EQUAL(p < r.locus(), true);
        }
        p = r.locus();
        first = false;
    }
}
