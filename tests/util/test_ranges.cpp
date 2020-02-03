#include "varoom/util/ranges.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ranges tests
#include <boost/test/unit_test.hpp>

#include <random>
#include <nlohmann/json.hpp>
#include "varoom/util/bitstream.hpp"
#include "varoom/util/elias_fano.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/json_builder.hpp"

namespace // anonymous
{
    uint64_t interlace(uint32_t x, uint32_t y)
    {
        uint64_t z = 0;
        uint64_t xb = 2;
        uint64_t yb = 1;
        while (x > 0 || y > 0)
        {
            if (x & 1)
            {
                z |= xb;
            }
            if (y & 1)
            {
                z |= yb;
            }

            x >>= 1;
            y >>= 1;
            xb <<= 2;
            yb <<= 2;
        }
        return z;
    }

    std::pair<uint32_t,uint32_t> deinterlace(uint64_t z)
    {
        uint32_t x = 0;
        uint32_t y = 0;
        uint32_t b = 1;
        while (z > 0)
        {
            if (z & 1)
            {
                y |= b;
            }
            z >>= 1;
            if (z & 1)
            {
                x |= b;
            }
            z >>= 1;
            b <<= 1;
        }
        return std::make_pair(x, y);
    }

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

    void build_random(uint64_t S, uint64_t D, uint64_t M, uint64_t N,
                      varoom::ranges_builder& R)
    {
        std::mt19937_64 rng(S);
        std::uniform_int_distribution<uint64_t> Z(0, D-1);
        std::poisson_distribution<uint64_t> W(M);
        std::vector<varoom::ranges_builder::range> R0;
        R0.reserve(N);
        for (uint64_t i = 0; i < N; ++i)
        {
            uint64_t x = Z(rng);
            uint64_t y = x + W(rng);
            while (y >= D)
            {
                x = Z(rng);
                y = x + W(rng);
            }
            R0.push_back(varoom::ranges_builder::range(x, y));
        }
        std::sort(R0.begin(), R0.end());
        for (uint64_t i = 0; i < R0.size(); ++i)
        {
            R.insert(R0[i]);
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
    R.insert(varoom::ranges_builder::range(5, 10));
    R.insert(varoom::ranges_builder::range(2, 6));

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

BOOST_AUTO_TEST_CASE( ranges_at_1 )
{
    varoom::ranges_builder R;
    R.insert(varoom::ranges_builder::range(5, 10));
    R.insert(varoom::ranges_builder::range(2, 6));

    std::vector<varoom::ranges_builder::range> rs;
    {
        rs.clear();
        R.ranges_at(0, 11, rs);
        BOOST_REQUIRE_EQUAL(rs.size(), 2);
        BOOST_CHECK_EQUAL(rs[0].first, 2);
        BOOST_CHECK_EQUAL(rs[0].second, 6);
        BOOST_CHECK_EQUAL(rs[1].first, 5);
        BOOST_CHECK_EQUAL(rs[1].second, 10);
    }
}

BOOST_AUTO_TEST_CASE( make_1 )
{
    varoom::ranges_builder R;
    R.insert(varoom::ranges_builder::range(5, 10));
    R.insert(varoom::ranges_builder::range(2, 6));

    varoom::ranges r;
    r.make(R);

    std::vector<varoom::ranges_builder::range> rs;
    {
        rs.clear();
        r.ranges_at(0, 11, rs);
        BOOST_REQUIRE_EQUAL(rs.size(), 2);
        BOOST_CHECK_EQUAL(rs[0].first, 2);
        BOOST_CHECK_EQUAL(rs[0].second, 6);
        BOOST_CHECK_EQUAL(rs[1].first, 5);
        BOOST_CHECK_EQUAL(rs[1].second, 10);
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
        std::vector<varoom::ranges::range> xs;
        R.ranges_at(bb, ee, xs);
        std::vector<varoom::ranges::range> ys;
        r.ranges_at(bb, ee, ys);

        BOOST_REQUIRE_EQUAL(ys.size(), xs.size());
        for (size_t j = 0; j < xs.size(); ++j)
        {
            BOOST_CHECK_EQUAL(ys[j].first, xs[j].first);
            BOOST_CHECK_EQUAL(ys[j].second, xs[j].second);
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
        std::vector<varoom::ranges::range> xs;
        R.ranges_at(bb, ee, xs);
        std::vector<varoom::ranges::range> ys;
        r.ranges_at(bb, ee, ys);

        BOOST_REQUIRE_EQUAL(ys.size(), xs.size());
        for (size_t j = 0; j < xs.size(); ++j)
        {
            BOOST_CHECK_EQUAL(ys[j].first, xs[j].first);
            BOOST_CHECK_EQUAL(ys[j].second, xs[j].second);
        }
    }
}

#if 0
BOOST_AUTO_TEST_CASE( succinct_ranges_2 )
{
    for (uint32_t n = 20; n <= 20; ++n)
    {
        uint64_t S = 17;
        uint64_t D = 1ULL << 32;
        uint64_t M = 1ULL << 12;
        uint64_t N = 1ULL << n;
        varoom::ranges_builder R;
        build_random(S, D, M, N, R);
    }

}
#endif
