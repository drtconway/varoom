#include "varoom/util/parsing.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE parsing tests
#include <boost/test/unit_test.hpp>

#include <boost/lexical_cast.hpp>

namespace // anonymous
{
    int wibble;

    std::unique_ptr<int> seventeen()
    {
        wibble += 1;
        return std::unique_ptr<int>(new int(17));
    }

    varoom::parsing::detail::stream<int> from(int p_first)
    {
        using cell_type = varoom::parsing::detail::cell<int>;
        using ptr_type = std::unique_ptr<cell_type>;
        using func_type = std::function<ptr_type()>;
        func_type f = [p_first]() {
            return ptr_type(new cell_type(p_first, from(p_first + 1)));
        };
        return varoom::parsing::detail::stream<int>(f);
    }

    varoom::parsing::detail::stream<int> between(int p_first, int p_last)
    {
        if (p_first > p_last)
        {
            return varoom::parsing::detail::stream<int>();
        }
        using cell_type = varoom::parsing::detail::cell<int>;
        using ptr_type = std::unique_ptr<cell_type>;
        using func_type = std::function<ptr_type()>;
        func_type f = [p_first, p_last]() {
            return ptr_type(new cell_type(p_first, between(p_first + 1, p_last)));
        };
        return varoom::parsing::detail::stream<int>(f);
    }

    struct simple_input : std::pair<std::string::const_iterator,std::string::const_iterator> {};

    typedef varoom::parsing::parse_result<simple_input,bool> simple_result;

    simple_result epsilon(simple_input p_in)
    {
        return varoom::parsing::detail::mreturn(std::pair<bool,simple_input>(true, p_in));
    }
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( sus1 )
{
    wibble = 1;
    varoom::parsing::detail::suspension<int> s(seventeen);
    BOOST_CHECK_EQUAL(wibble, 1);
    BOOST_CHECK_EQUAL(s.get(), 17);
    BOOST_CHECK_EQUAL(wibble, 2);
    BOOST_CHECK_EQUAL(s.get(), 17);
    BOOST_CHECK_EQUAL(wibble, 2);
}

BOOST_AUTO_TEST_CASE( stream1 )
{
    varoom::parsing::detail::stream<int> s = from(1).take(10);

    std::vector<int> xs;
    varoom::parsing::detail::for_each(s, [&](int x) {
        xs.push_back(x);
    });

    BOOST_REQUIRE_EQUAL(xs.size(), 10);
    BOOST_CHECK_EQUAL(xs[0], 1);
    BOOST_CHECK_EQUAL(xs[1], 2);
    BOOST_CHECK_EQUAL(xs[8], 9);
    BOOST_CHECK_EQUAL(xs[9], 10);
}

BOOST_AUTO_TEST_CASE( stream2 )
{
    varoom::parsing::detail::stream<int> s = between(1, 10);

    std::vector<int> xs;
    varoom::parsing::detail::for_each(s, [&](int x) {
        xs.push_back(x);
    });

    BOOST_REQUIRE_EQUAL(xs.size(), 10);
    BOOST_CHECK_EQUAL(xs[0], 1);
    BOOST_CHECK_EQUAL(xs[1], 2);
    BOOST_CHECK_EQUAL(xs[8], 9);
    BOOST_CHECK_EQUAL(xs[9], 10);
}

BOOST_AUTO_TEST_CASE( fmap1 )
{
    varoom::parsing::detail::stream<int> s = from(1).take(10);

    std::function<std::string(int)> f = [](int x) { return boost::lexical_cast<std::string>(x); };

    varoom::parsing::detail::stream<std::string> t = varoom::parsing::detail::fmap(s, f);

    std::vector<std::string> xs;
    varoom::parsing::detail::for_each(t, [&](const std::string& x) {
        xs.push_back(x);
    });

    BOOST_REQUIRE_EQUAL(xs.size(), 10);
    BOOST_CHECK_EQUAL(xs[0], "1");
    BOOST_CHECK_EQUAL(xs[1], "2");
    BOOST_CHECK_EQUAL(xs[8], "9");
    BOOST_CHECK_EQUAL(xs[9], "10");
}

BOOST_AUTO_TEST_CASE( mbind1 )
{
    varoom::parsing::detail::stream<int> s = varoom::parsing::detail::mbind(from(1), [](int x) { return between(1, x); }).take(10);

    std::vector<int> xs;
    varoom::parsing::detail::for_each(s, [&](const int& x) {
        xs.push_back(x);
    });

    BOOST_REQUIRE_EQUAL(xs.size(), 10);
    BOOST_CHECK_EQUAL(xs[0], 1);
    BOOST_CHECK_EQUAL(xs[1], 1);
    BOOST_CHECK_EQUAL(xs[8], 3);
    BOOST_CHECK_EQUAL(xs[9], 4);
}

BOOST_AUTO_TEST_CASE( parser1 )
{
    varoom::parsing::parser<simple_input,bool> p = epsilon;
}
