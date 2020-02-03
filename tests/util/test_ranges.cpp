#include "varoom/util/ranges.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ranges tests
#include <boost/test/unit_test.hpp>

#include <random>
#include <boost/lexical_cast.hpp>
#include "varoom/util/bitstream.hpp"
#include "varoom/util/elias_fano.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/json_builder.hpp"

namespace // anonymous
{
#if 0
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
#endif

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

#if 0
    uint64_t build_random(uint64_t S, uint64_t D, uint64_t M, uint64_t N,
                      varoom::ranges_builder& R)
    {
        std::mt19937_64 rng(S);
        std::uniform_int_distribution<uint64_t> Z(0, D-1);
        //std::poisson_distribution<uint64_t> W(M);
        std::exponential_distribution<> W(1.0);
        std::vector<varoom::ranges_builder::range> R0;
        R0.reserve(N);
        uint64_t t = 0;
        for (uint64_t i = 0; i < N; ++i)
        {
            uint64_t x = Z(rng);
            //uint64_t y = x + M;
            //uint64_t y = x + W(rng);
            uint64_t y = x + M * (1 + W(rng));
            while (y >= D)
            {
                x = Z(rng);
                //y = x + M;
                //y = x + W(rng);
                y = x + M * (1 + W(rng));
            }
            t += y - x;
            R0.push_back(varoom::ranges_builder::range(x, y));
        }
        std::sort(R0.begin(), R0.end());
        for (uint64_t i = 0; i < R0.size(); ++i)
        {
            R.insert(R0[i]);
        }
        return t;
    }
#endif
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
    //
    // Not really a test case
    //

    size_t d = 17;
    for (uint32_t i = 0; i < 100; ++i)
    {
        for (uint32_t n = 16; n < 22; ++n, ++d)
        {
            uint64_t S = d;
            uint64_t D = 1ULL << 32;
            uint64_t M = 1ULL << 10;
            uint64_t N = 1ULL << n;
            varoom::ranges_builder R;
            uint64_t t = build_random(S, D, M, N, R);

            varoom::ranges r;
            r.make(R);

            nlohmann::json s = r.stats();
            s["C"] = t;
            s["S"] = S;
            s["N"] = N;

            size_t mB = s["B"]["memory"];
            size_t mE = s["E"]["memory"];
            size_t mTi = s["Ti"]["memory"];
            size_t mTe = s["Te"]["memory"];
            size_t m = mB + mE + mTi + mTe;

            std::vector<uint64_t> h = s["Ti"]["hist"];

            uint64_t hs = 0;
            uint64_t hn = 0;
            uint64_t hm = h.size() - 1;
            for (size_t j = 0; j < h.size(); ++j)
            {
                if (h[j] == 0)
                {
                    continue;
                }
                hs += j*h[j];
                hn += h[j];
            }

            std::cerr << S
                << '\t' << n
                << '\t' << N
                << '\t' << (double(t)/double(D))
                << '\t' << (double(s["Ti"]["count"])/double(s["Ti"]["size"]))
                << '\t' << (double(hs)/double(hn))
                << '\t' << hm
                << '\t' << (double(m)/double(1024*1024))
                << '\t' << (double(m)/double(N))
                << std::endl;
        }
    }
    return;
    for (uint32_t n = 16; n < 22; ++n)
    {
        uint64_t S = 17 + n;
        uint64_t D = 1ULL << 32;
        uint64_t M = 1ULL << 12;
        uint64_t N = 1ULL << n;
        varoom::ranges_builder R;
        build_random(S, D, M, N, R);

        varoom::ranges r;
        r.make(R);

        nlohmann::json s = r.stats();
        s["N"] = N;
        std::cerr << s << std::endl;

        std::string fn = std::string("tests/tmp/rnd-") + boost::lexical_cast<std::string>(N) + std::string(".sdsl");
        varoom::output_file_holder_ptr outp = varoom::files::out(fn);
        r.save(**outp);
    }
}
#endif
