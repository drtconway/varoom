#include "varoom/util/gdbm.hpp"
#include "varoom/util/files.hpp"
#include <boost/lexical_cast.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE gdbm tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
    struct foo
    {
        int baz;
        double qux;
    };
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( try1 )
{
    std::string nm = "tests/tmp/foo.gdbm";
    varoom::files::remove(nm);
    {
        varoom::gdbm G(nm, varoom::gdbm::trunc);
        BOOST_CHECK_EQUAL(G.size(), 0);
        BOOST_CHECK_EQUAL(G.count("bar"), false);
        G.put("bar", "baz");
        BOOST_CHECK_EQUAL(G.size(), 1);
        BOOST_CHECK_EQUAL(G.count("bar"), true);
    }
    {
        varoom::gdbm G(nm, varoom::gdbm::read_only);
        BOOST_CHECK_EQUAL(G.size(), 1);
        BOOST_CHECK_EQUAL(G.count("bar"), true);
        BOOST_CHECK_EQUAL(G.get("bar").str(), "baz");
    }
}

BOOST_AUTO_TEST_CASE( try2 )
{
    std::string nm = "tests/tmp/foo.gdbm";
    const size_t N = 1000;
    varoom::files::remove(nm);
    std::map<std::string,std::string> m;
    {
        varoom::gdbm G(nm, varoom::gdbm::trunc);
        for (size_t i = 0; i < N; ++i)
        {
            double j = std::sin(6.28 * double(i) / double(N));
            std::string k = boost::lexical_cast<std::string>(i);
            std::string v = boost::lexical_cast<std::string>(j);
            G.put(k, v);
            m[k] = v;
        }
    }
    {
        varoom::gdbm G(nm, varoom::gdbm::read_only);
        BOOST_CHECK_EQUAL(G.size(), N);
        varoom::gdbm::owning_datum x;
        while (G.next(x))
        {
            std::string k = x.str();
            BOOST_CHECK_EQUAL(m.count(k), true);
            BOOST_CHECK_EQUAL(G.count(k), true);
            m.erase(k);
        }
        BOOST_CHECK_EQUAL(m.size(), 0);
    }
}

BOOST_AUTO_TEST_CASE( try3 )
{
    std::string nm = "tests/tmp/foo.gdbm";
    varoom::files::remove(nm);
    {
        varoom::gdbm G(nm, varoom::gdbm::trunc);
        std::string k = "bar";
        foo v = { -42, 3.1415 };
        G.put(varoom::gdbm::str_datum(k), varoom::gdbm::make_datum(v));
    }
    {
        varoom::gdbm G(nm, varoom::gdbm::read_only);
        std::string k = "bar";
        foo v = varoom::gdbm::datum_cast<foo>(G.get(k));
        BOOST_CHECK_EQUAL(v.baz, -42);
        BOOST_CHECK_EQUAL(v.qux, 3.1415);
    }
}
