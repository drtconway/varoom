#include "varoom/seq/compact_seq.hpp"

#include <nlohmann/json.hpp>
#include "varoom/util/files.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE compact_seq tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( test_create )
{
    if (0)
    {
        varoom::seq::compact_seq::make("data/chr13.fa.gz", "data/wibble");

        varoom::seq::compact_seq chr13("data/wibble");

        std::cerr << "getting..." << std::endl;
        std::string s = chr13.get("chr13", 32890540, 32890720);
        std::cout << s << std::endl;
    }
}

