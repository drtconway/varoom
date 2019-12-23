#include "varoom/seq/compact.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE codec8 tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( enc1 )
{
    varoom::seq::compact::make_csa("data/hg19/chr1.txt", "data/compact/chr1");
}
