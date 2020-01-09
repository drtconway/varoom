#include "varoom/util/ranges.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ranges tests
#include <boost/test/unit_test.hpp>

#include <nlohmann/json.hpp>
#include "varoom/util/files.hpp"
#include "varoom/util/json_builder.hpp"

namespace // anonymous
{
    void build_chr1(varoom::ranges_builder& R)
    {
        std::ifstream i("tests/data/chr1-exons.json");
        nlohmann::json j;
        i >> j;
        for (size_t k = 0; k < j.size(); ++k)
        {
            uint64_t p0 = j[k][0];
            uint64_t p1 = j[k][1];
            R.insert(varoom::ranges_builder::range(p0, p1));
        }
    }
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( simple1 )
{
    varoom::ranges_builder R;
    R.insert(varoom::ranges_builder::range(6, 7));
    R.insert(varoom::ranges_builder::range(1, 4));
    R.insert(varoom::ranges_builder::range(8, 12));

    std::vector<varoom::ranges_builder::range> rs;
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
    varoom::ranges_builder R;
    auto r0 = R.insert(varoom::ranges_builder::range(5, 10));
    auto r1 = R.insert(varoom::ranges_builder::range(2, 6));

    BOOST_CHECK_EQUAL(r0, 0);
    BOOST_CHECK_EQUAL(r1, 1);

    std::vector<varoom::ranges_builder::range> rs;
    {
        rs.clear();
        R.overlapping_ranges(0, 10, rs);
        BOOST_REQUIRE_EQUAL(rs.size(), 3);
        BOOST_CHECK_EQUAL(rs[0].first, 2);
        BOOST_CHECK_EQUAL(rs[0].second, 5);
        BOOST_REQUIRE_EQUAL(R.ranges_at_segment(2).size(), 1);
        BOOST_CHECK_EQUAL(R.ranges_at_segment(2)[0], 1);
        BOOST_CHECK_EQUAL(rs[1].first, 5);
        BOOST_CHECK_EQUAL(rs[1].second, 6);
        BOOST_REQUIRE_EQUAL(R.ranges_at_segment(5).size(), 2);
        BOOST_CHECK_EQUAL(R.ranges_at_segment(5)[0], 0);
        BOOST_CHECK_EQUAL(R.ranges_at_segment(5)[1], 1);
        BOOST_CHECK_EQUAL(rs[2].first, 6);
        BOOST_CHECK_EQUAL(rs[2].second, 10);
        BOOST_REQUIRE_EQUAL(R.ranges_at_segment(6).size(), 1);
        BOOST_CHECK_EQUAL(R.ranges_at_segment(6)[0], 0);
    }
}

BOOST_AUTO_TEST_CASE( succinct_ranges_1 )
{
    varoom::ranges_builder R;
    build_chr1(R);

    if (0)
    {
        const std::unordered_map<size_t,std::vector<uint32_t>>& idx = R.index();

        varoom::json_writer J(std::cout);
        J.begin_array();
        for (auto itr = idx.begin(); itr != idx.end(); ++itr)
        {
            J.begin_array();
            J.number_unsigned(itr->first);
            J.begin_array();
            for (size_t i = 0; i < itr->second.size(); ++i)
            {
                J.number_unsigned(itr->second[i]);
            }
            J.end_array();
            J.end_array();
        }
        J.end_array();
    }

    varoom::ranges r;
    r.make(R);

    size_t N = R.begins().size();
    for (size_t i = 0; i < N; ++i)
    {
        size_t bb = R.begins().select(i);
        size_t ee = R.ends().select(i);
        std::vector<varoom::ranges::range_id> xs;
        R.ranges_at(bb, ee, xs);
        std::vector<varoom::ranges::range_id> ys;
        r.ranges_at(bb, ee, ys);

        BOOST_REQUIRE_EQUAL(ys.size(), xs.size());
        for (size_t j = 0; j < xs.size(); ++j)
        {
            BOOST_CHECK_EQUAL(ys[j], xs[j]);
        }
    }

    {
        varoom::output_file_holder_ptr outp = varoom::files::out("tests/tmp/foo.gz");
        r.save(**outp);
    }
    {
        varoom::input_file_holder_ptr inp = varoom::files::in("tests/tmp/foo.gz");
        r.load(**inp);
    }

    for (size_t i = 0; i < N; ++i)
    {
        size_t bb = R.begins().select(i);
        size_t ee = R.ends().select(i);
        std::vector<varoom::ranges::range_id> xs;
        R.ranges_at(bb, ee, xs);
        std::vector<varoom::ranges::range_id> ys;
        r.ranges_at(bb, ee, ys);

        BOOST_REQUIRE_EQUAL(ys.size(), xs.size());
        for (size_t j = 0; j < xs.size(); ++j)
        {
            BOOST_CHECK_EQUAL(ys[j], xs[j]);
        }
    }
}
