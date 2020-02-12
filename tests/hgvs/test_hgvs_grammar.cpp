#include "varoom/hgvs/hgvs_grammar.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE codec8 tests
#include <boost/test/unit_test.hpp>

using namespace varoom::funds;

namespace // anonymous
{
    template <typename T>
    size_t length(list<T> xs)
    {
        size_t n = 0;
        while (!xs.empty())
        {
            xs = xs.tail();
            n += 1;
        }
        return n;
    }

    template <typename T>
    T access(list<T> xs, size_t i)
    {
        while (i > 0)
        {
            xs = xs.tail();
            i -= 1;
        }
        return xs.head();
    }
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( digits1 )
{
    using P = parser<uint64_t>;
    P p = varoom::hgvs::grammar::num();
    {
        auto s = parse(p, symbols("2"));
        BOOST_REQUIRE_EQUAL(length(s), 1);
        BOOST_CHECK_EQUAL(access(s, 0).first, 2);
    }
}

BOOST_AUTO_TEST_CASE( acc1 )
{
    using P = parser<varoom::hgvs::accession>;
    P p = varoom::hgvs::grammar::accn();
    auto s = parse(p, symbols("chr1"));
    BOOST_REQUIRE_EQUAL(length(s), 1);
    BOOST_CHECK_EQUAL(access(s, 0).first.name, "chr1");
    BOOST_CHECK_EQUAL(access(s, 0).first.version, 0);
}

BOOST_AUTO_TEST_CASE( acc2 )
{
    using P = parser<varoom::hgvs::accession>;
    P p = varoom::hgvs::grammar::accn();
    auto s = parse(p, symbols("NC_001802"));
    BOOST_REQUIRE_EQUAL(length(s), 1);
    BOOST_CHECK_EQUAL(access(s, 0).first.name, "NC_001802");
    BOOST_CHECK_EQUAL(access(s, 0).first.version, 0);
}

BOOST_AUTO_TEST_CASE( acc3 )
{
    using P = parser<varoom::hgvs::accession>;
    P p = varoom::hgvs::grammar::accn();
    auto s = parse(p, symbols("NC_001802.1"));
    BOOST_REQUIRE_EQUAL(length(s), 1);
    BOOST_CHECK_EQUAL(access(s, 0).first.name, "NC_001802");
    BOOST_CHECK_EQUAL(access(s, 0).first.version, 1);
}

BOOST_AUTO_TEST_CASE( hgvsg_id_1 )
{
    using P = parser<varoom::hgvs::variant_ptr>;
    P p = varoom::hgvs::grammar::g_var();
    auto s = parse(p, symbols("NC_000023.10:g.33038255="));
    BOOST_REQUIRE_EQUAL(length(s), 1);
    varoom::hgvs::variant_ptr v = access(s, 0).first;
    const varoom::hgvs::hgvsg_id* w = dynamic_cast<const varoom::hgvs::hgvsg_id*>(v.get());
    BOOST_REQUIRE_EQUAL(w != NULL, true);
    BOOST_CHECK_EQUAL(w->acc.name, "NC_000023");
    BOOST_CHECK_EQUAL(w->acc.version, 10);
    BOOST_CHECK_EQUAL(static_cast<uint64_t>(w->pos), 33038255);
}

BOOST_AUTO_TEST_CASE( hgvsg_sub_1 )
{
    using P = parser<varoom::hgvs::variant_ptr>;
    P p = varoom::hgvs::grammar::g_var();
    auto s = parse(p, symbols("NC_000023.10:g.33038255C>A"));
    BOOST_REQUIRE_EQUAL(length(s), 1);
    varoom::hgvs::variant_ptr v = access(s, 0).first;
    const varoom::hgvs::hgvsg_sub* w = dynamic_cast<const varoom::hgvs::hgvsg_sub*>(v.get());
    BOOST_REQUIRE_EQUAL(w != NULL, true);
    BOOST_CHECK_EQUAL(w->acc.name, "NC_000023");
    BOOST_CHECK_EQUAL(w->acc.version, 10);
    BOOST_CHECK_EQUAL(static_cast<uint64_t>(w->pos), 33038255);
    BOOST_CHECK_EQUAL(w->ref, 'C');
    BOOST_CHECK_EQUAL(w->alt, 'A');
}

BOOST_AUTO_TEST_CASE( hgvsg_del_1 )
{
    using P = parser<varoom::hgvs::variant_ptr>;
    P p = varoom::hgvs::grammar::g_var();
    auto s = parse(p, symbols("NG_012232.1:g.19del"));
    BOOST_REQUIRE_EQUAL(length(s), 1);
    varoom::hgvs::variant_ptr v = access(s, 0).first;
    const varoom::hgvs::hgvsg_del* w = dynamic_cast<const varoom::hgvs::hgvsg_del*>(v.get());
    BOOST_REQUIRE_EQUAL(w != NULL, true);
    BOOST_CHECK_EQUAL(w->acc.name, "NG_012232");
    BOOST_CHECK_EQUAL(w->acc.version, 1);
    BOOST_CHECK_EQUAL(static_cast<uint64_t>(w->loc.first), 19);
    BOOST_CHECK_EQUAL(static_cast<uint64_t>(w->loc.second), 19);
    BOOST_CHECK_EQUAL(w->ref, "");
}

BOOST_AUTO_TEST_CASE( hgvsg_del_2 )
{
    using P = parser<varoom::hgvs::variant_ptr>;
    P p = varoom::hgvs::grammar::g_var();
    auto s = parse(p, symbols("NG_012232.1:g.19delT"));
    BOOST_REQUIRE_EQUAL(length(s), 1);
    varoom::hgvs::variant_ptr v = access(s, 0).first;
    const varoom::hgvs::hgvsg_del* w = dynamic_cast<const varoom::hgvs::hgvsg_del*>(v.get());
    BOOST_REQUIRE_EQUAL(w != NULL, true);
    BOOST_CHECK_EQUAL(w->acc.name, "NG_012232");
    BOOST_CHECK_EQUAL(w->acc.version, 1);
    BOOST_CHECK_EQUAL(static_cast<uint64_t>(w->loc.first), 19);
    BOOST_CHECK_EQUAL(static_cast<uint64_t>(w->loc.second), 19);
    BOOST_CHECK_EQUAL(w->ref, "T");
}

BOOST_AUTO_TEST_CASE( hgvsg_del_3 )
{
    using P = parser<varoom::hgvs::variant_ptr>;
    P p = varoom::hgvs::grammar::g_var();
    auto s = parse(p, symbols("NG_012232.1:g.19_21del"));
    BOOST_REQUIRE_EQUAL(length(s), 1);
    varoom::hgvs::variant_ptr v = access(s, 0).first;
    const varoom::hgvs::hgvsg_del* w = dynamic_cast<const varoom::hgvs::hgvsg_del*>(v.get());
    BOOST_REQUIRE_EQUAL(w != NULL, true);
    BOOST_CHECK_EQUAL(w->acc.name, "NG_012232");
    BOOST_CHECK_EQUAL(w->acc.version, 1);
    BOOST_CHECK_EQUAL(static_cast<uint64_t>(w->loc.first), 19);
    BOOST_CHECK_EQUAL(static_cast<uint64_t>(w->loc.second), 21);
    BOOST_CHECK_EQUAL(w->ref, "");
}

BOOST_AUTO_TEST_CASE( hgvsg_del_4 )
{
    using P = parser<varoom::hgvs::variant_ptr>;
    P p = varoom::hgvs::grammar::g_var();
    auto s = parse(p, symbols("NG_012232.1:g.19_21delTCA"));
    BOOST_REQUIRE_EQUAL(length(s), 1);
    varoom::hgvs::variant_ptr v = access(s, 0).first;
    const varoom::hgvs::hgvsg_del* w = dynamic_cast<const varoom::hgvs::hgvsg_del*>(v.get());
    BOOST_REQUIRE_EQUAL(w != NULL, true);
    BOOST_CHECK_EQUAL(w->acc.name, "NG_012232");
    BOOST_CHECK_EQUAL(w->acc.version, 1);
    BOOST_CHECK_EQUAL(static_cast<uint64_t>(w->loc.first), 19);
    BOOST_CHECK_EQUAL(static_cast<uint64_t>(w->loc.second), 21);
    BOOST_CHECK_EQUAL(w->ref, "TCA");
}
