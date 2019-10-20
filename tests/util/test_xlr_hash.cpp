#include "varoom/util/xlr_hash.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE xlr_hash tests
#include <boost/test/unit_test.hpp>

namespace
{
    class bit_counter
    {
    public:
        void add(std::uint64_t p_x)
        {
        }
    };
}

BOOST_AUTO_TEST_CASE( basic )
{
    using namespace std;
    using namespace varoom;
    
    const uint64_t s = 0;
    const uint64_t x = 0;
    BOOST_CHECK_EQUAL(xlr_hash::hash(s, x), 0xa79f316e67b0852bULL);
}

BOOST_AUTO_TEST_CASE( basic )
{
    using namespace std;
    using namespace varoom;
    
    const uint64_t s = 0;
    const uint64_t x = 0;
    BOOST_CHECK_EQUAL(xlr_hash::hash(s, x), 0xa79f316e67b0852bULL);
}
