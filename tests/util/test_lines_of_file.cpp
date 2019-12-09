#include "varoom/util/lines_of_file.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE lines_of_file tests
#include <boost/test/unit_test.hpp>

#include <iostream>
#include "varoom/util/files.hpp"

namespace // anonymous
{
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( lines1 )
{
    varoom::input_file_holder_ptr inp = varoom::files::in("tests/data/SRR2912526-100_1.fastq.gz");
    
    varoom::lines_of_file<> l(**inp);

    std::string s;
    int n = 0;
    while (l.next(s))
    {
        ++n;
    }
    BOOST_CHECK_EQUAL(n, 100);
}

BOOST_AUTO_TEST_CASE( lines2 )
{
    varoom::input_file_holder_ptr inp = varoom::files::in("tests/data/SRR2912526-100_1.fastq.gz");
    
    varoom::lines_of_file<false> l(**inp);

    std::string s;
    int n = 0;
    while (l.next(s))
    {
        ++n;
    }
    BOOST_CHECK_EQUAL(n, 100);
}

BOOST_AUTO_TEST_CASE( lines3 )
{
    varoom::input_file_holder_ptr inp = varoom::files::in("tests/data/SRR2912526-100_1.fastq.gz");
    
    varoom::lines_of_file<true> l(**inp);

    std::string s;
    int n = 0;
    while (l.next(s))
    {
        ++n;
    }
    BOOST_CHECK_EQUAL(n, 100);
}
