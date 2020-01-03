#include "varoom/util/ranges.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ranges tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( simple1 )
{
    varoom::ranges R;
    R.insert(varoom::ranges::range(6, 7));
    R.insert(varoom::ranges::range(1, 4));
    R.insert(varoom::ranges::range(8, 12));

    std::vector<varoom::ranges::range> rs;
    {
        rs.clear();
        R.overlapping_ranges(0, 10, rs);
        BOOST_REQUIRE_EQUAL(rs.size(), 3);
        BOOST_CHECK_EQUAL(rs[0].first, 1);
        BOOST_CHECK_EQUAL(rs[0].second, 4);
        BOOST_CHECK_EQUAL(rs[1].first, 6);
        BOOST_CHECK_EQUAL(rs[1].second, 7);
        BOOST_CHECK_EQUAL(rs[2].first, 8);
        BOOST_CHECK_EQUAL(rs[2].second, 12);
    }
    {
        rs.clear();
        R.overlapping_ranges(2, 9, rs);
        BOOST_REQUIRE_EQUAL(rs.size(), 3);
        BOOST_CHECK_EQUAL(rs[0].first, 1);
        BOOST_CHECK_EQUAL(rs[0].second, 4);
        BOOST_CHECK_EQUAL(rs[1].first, 6);
        BOOST_CHECK_EQUAL(rs[1].second, 7);
        BOOST_CHECK_EQUAL(rs[2].first, 8);
        BOOST_CHECK_EQUAL(rs[2].second, 12);
    }
    {
        rs.clear();
        R.overlapping_ranges(4, 6, rs);
        BOOST_REQUIRE_EQUAL(rs.size(), 0);
    }
    {
        rs.clear();
        R.overlapping_ranges(5, 8, rs);
        BOOST_REQUIRE_EQUAL(rs.size(), 1);
        BOOST_CHECK_EQUAL(rs[0].first, 6);
        BOOST_CHECK_EQUAL(rs[0].second, 7);
    }
    {
        rs.clear();
        R.overlapping_ranges(1, 4, rs);
        BOOST_REQUIRE_EQUAL(rs.size(), 1);
        BOOST_CHECK_EQUAL(rs[0].first, 1);
        BOOST_CHECK_EQUAL(rs[0].second, 4);
    }
    {
        rs.clear();
        R.overlapping_ranges(6, 7, rs);
        BOOST_REQUIRE_EQUAL(rs.size(), 1);
        BOOST_CHECK_EQUAL(rs[0].first, 6);
        BOOST_CHECK_EQUAL(rs[0].second, 7);
    }
    {
        rs.clear();
        R.overlapping_ranges(8, 12, rs);
        BOOST_REQUIRE_EQUAL(rs.size(), 1);
        BOOST_CHECK_EQUAL(rs[0].first, 8);
        BOOST_CHECK_EQUAL(rs[0].second, 12);
    }
}

BOOST_AUTO_TEST_CASE( split1 )
{
    varoom::ranges R;
    auto r0 = R.insert(varoom::ranges::range(5, 10));
    auto r1 = R.insert(varoom::ranges::range(2, 6));

    BOOST_CHECK_EQUAL(r0, 0);
    BOOST_CHECK_EQUAL(r1, 1);

    std::vector<varoom::ranges::range> rs;
    {
        rs.clear();
        R.overlapping_ranges(0, 10, rs);
        BOOST_REQUIRE_EQUAL(rs.size(), 3);
        BOOST_CHECK_EQUAL(rs[0].first, 2);
        BOOST_CHECK_EQUAL(rs[0].second, 5);
        BOOST_REQUIRE_EQUAL(R.ranges_at(2).size(), 1);
        BOOST_CHECK_EQUAL(R.ranges_at(2)[0], 1);
        BOOST_CHECK_EQUAL(rs[1].first, 5);
        BOOST_CHECK_EQUAL(rs[1].second, 6);
        BOOST_REQUIRE_EQUAL(R.ranges_at(5).size(), 2);
        BOOST_CHECK_EQUAL(R.ranges_at(5)[0], 0);
        BOOST_CHECK_EQUAL(R.ranges_at(5)[1], 1);
        BOOST_CHECK_EQUAL(rs[2].first, 6);
        BOOST_CHECK_EQUAL(rs[2].second, 10);
        BOOST_REQUIRE_EQUAL(R.ranges_at(6).size(), 1);
        BOOST_CHECK_EQUAL(R.ranges_at(6)[0], 0);
    }
}
