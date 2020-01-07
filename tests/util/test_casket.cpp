#include "varoom/util/casket.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE casket tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
    void check_uint64_at(const std::string& p_str, size_t p_pos, uint64_t x)
    {
        BOOST_CHECK_EQUAL(p_str.size() >= p_pos + sizeof(x), true);
        const char* p = &p_str[p_pos];
        const char* q = reinterpret_cast<const char*>(&x);
        for (size_t i = 0; i < sizeof(x); ++i, ++p, ++q)
        {
            BOOST_CHECK_EQUAL(*p, *q);
        }
    }
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( ranged_1 )
{
    std::string s = "0123456789abcdef";
    varoom::detail::bytes_input_file_holder i(s);
    varoom::detail::ranged_input_file_holder r(*i, 8, 4);
    std::string x("....");
    (*r).read(&x[0], 4);
    BOOST_CHECK_EQUAL(x, "89ab");
}

BOOST_AUTO_TEST_CASE( bytes_out_tellp_1 )
{
    std::string s = "0123456789abcdef";
    std::string t;
    std::function<void(const std::vector<char>&)> f = [&](const std::vector<char>& v) mutable {
        t.clear();
        t.insert(t.end(), v.begin(), v.end());
    };
    {
        varoom::detail::bytes_output_file_holder r(f);
        // BOOST_CHECK_EQUAL((*r).tellp(), 0);
        (*r).write(&s[0], s.size());
        // BOOST_CHECK_EQUAL((*r).tellp(), s.size());
    }
    BOOST_CHECK_EQUAL(t, s);
}

BOOST_AUTO_TEST_CASE( bytes_in_tellg_1 )
{
    std::string s = "0123456789abcdef";
    varoom::detail::bytes_input_file_holder i(s);
    BOOST_CHECK_EQUAL((*i).tellg(), 0);
    uint64_t x;
    (*i).read(reinterpret_cast<char*>(&x), sizeof(x));
    BOOST_CHECK_EQUAL((*i).tellg(), sizeof(x));
}

BOOST_AUTO_TEST_CASE( bytes_in_seekg_1 )
{
    std::string s = "0123456789abcdef";
    varoom::detail::bytes_input_file_holder i(s);
    BOOST_CHECK_EQUAL((*i).tellg(), 0);
    (*i).seekg(8, std::ios_base::beg);
    BOOST_CHECK_EQUAL((*i).tellg(), 8);
}

BOOST_AUTO_TEST_CASE( bytes_in_seekg_2 )
{
    std::string s = "0123456789abcdef";
    varoom::detail::bytes_input_file_holder i(s);
    BOOST_CHECK_EQUAL((*i).tellg(), 0);
    (*i).seekg(0, std::ios_base::end);
    uint64_t l = (*i).tellg();
    BOOST_CHECK_EQUAL(l, 16);
    (*i).seekg(l - 8, std::ios_base::beg);
    BOOST_CHECK_EQUAL((*i).tellg(), 8);
}

BOOST_AUTO_TEST_CASE( file_in_seekg_2 )
{
    std::string s = "0123456789abcdef";
    {
        std::ofstream o("tests/tmp/foo", std::ios_base::binary);
        o.write(&s[0], s.size());
    }
    std::ifstream i("tests/tmp/foo");
    BOOST_CHECK_EQUAL(i.tellg(), 0);
    i.seekg(0, std::ios_base::end);
    uint64_t l = i.tellg();
    BOOST_CHECK_EQUAL(l, 16);
    i.seekg(l - 8, std::ios_base::beg);
    BOOST_CHECK_EQUAL(i.tellg(), 8);
}

BOOST_AUTO_TEST_CASE( empty_casket_write )
{
    std::string t;
    std::function<void(const std::vector<char>&)> f = [&](const std::vector<char>& v) mutable {
        t.clear();
        t.insert(t.end(), v.begin(), v.end());
    };
    {
        varoom::detail::bytes_output_file_holder r(f);
        varoom::casket c(*r);
    }
    BOOST_CHECK_EQUAL(t.size(), 16);
    check_uint64_at(t, 0, 0);
    check_uint64_at(t, 8, 0);
}

BOOST_AUTO_TEST_CASE( empty_casket_read )
{
    std::string t;
    std::function<void(const std::vector<char>&)> f = [&](const std::vector<char>& v) mutable {
        t.clear();
        t.insert(t.end(), v.begin(), v.end());
    };
    {
        varoom::detail::bytes_output_file_holder r(f);
        varoom::casket c(*r);
    }
    BOOST_CHECK_EQUAL(t.size(), 16);
    check_uint64_at(t, 0, 0);
    check_uint64_at(t, 8, 0);

    varoom::detail::bytes_input_file_holder i(t);
    varoom::casket c(*i);
    std::vector<std::string> v = c.contents();
    BOOST_REQUIRE_EQUAL(v.size(), 0);
}

BOOST_AUTO_TEST_CASE( small_casket )
{
    std::string s("0123456789abcdef");
    std::string t;
    std::function<void(const std::vector<char>&)> f = [&](const std::vector<char>& v) mutable {
        t.clear();
        t.insert(t.end(), v.begin(), v.end());
    };
    {
        varoom::detail::bytes_output_file_holder r(f);
        varoom::casket c(*r);
        varoom::output_file_holder_ptr op = c.out("wibble-1");
        (**op).write(&s[0], s.size());
    }
    // Expected size:
    //   16     for the string in s;
    // + 8      for the toc len;
    // + 8      for the offset of wibble-1
    // + 8      for the length of wibble-1
    // + 8      for the length of the name "wibble-1"
    // + 8      for the characters "wibble-1"
    // + 8      for the pointer to the toc
    // = 64
    BOOST_REQUIRE_EQUAL(t.size(), 64);
    BOOST_CHECK_EQUAL(t.substr(0, 16), s);
    check_uint64_at(t, 16, 1);
    check_uint64_at(t, 24, 0);
    check_uint64_at(t, 32, 16);
    check_uint64_at(t, 40, 8);
    BOOST_CHECK_EQUAL(t.substr(48, 8), "wibble-1");
    check_uint64_at(t, 56, 16);

    varoom::detail::bytes_input_file_holder i(t);
    varoom::casket c(*i);
    std::vector<std::string> v = c.contents();
    BOOST_REQUIRE_EQUAL(v.size(), 1);
    BOOST_REQUIRE_EQUAL(v[0], "wibble-1");
    {
        varoom::input_file_holder_ptr ip = c.in("wibble-1");
        std::string u;
        u.resize(16);
        (**ip).read(&u[0], 16);
        BOOST_CHECK_EQUAL(u, s);
    }
}

BOOST_AUTO_TEST_CASE( medium_casket )
{
    std::string s0("0123456789abcdef");
    std::string s1("abcdef0123456789");
    std::string t;
    std::function<void(const std::vector<char>&)> f = [&](const std::vector<char>& v) mutable {
        t.clear();
        t.insert(t.end(), v.begin(), v.end());
    };
    {
        varoom::detail::bytes_output_file_holder r(f);
        varoom::casket c(*r);
        {
            varoom::output_file_holder_ptr op = c.out("wibble-1");
            (**op).write(&s0[0], s0.size());
        }
        {
            varoom::output_file_holder_ptr op = c.out("wibble-2");
            (**op).write(&s1[0], s1.size());
        }
    }
    // Expected size:
    //   16     for the string in s0;
    //   16     for the string in s1;
    // + 8      for the toc len;
    // + 8      for the offset of wibble-1
    // + 8      for the length of wibble-1
    // + 8      for the length of the name "wibble-1"
    // + 8      for the characters "wibble-1"
    // + 8      for the offset of wibble-2
    // + 8      for the length of wibble-2
    // + 8      for the length of the name "wibble-2"
    // + 8      for the characters "wibble-2"
    // + 8      for the pointer to the toc
    // = 64
    BOOST_REQUIRE_EQUAL(t.size(), 112);
    BOOST_CHECK_EQUAL(t.substr(0, 16), s0);
    BOOST_CHECK_EQUAL(t.substr(16, 16), s1);
    check_uint64_at(t, 32, 2);
    check_uint64_at(t, 40, 0);
    check_uint64_at(t, 48, 16);
    check_uint64_at(t, 56, 8);
    BOOST_CHECK_EQUAL(t.substr(64, 8), "wibble-1");
    check_uint64_at(t, 72, 16);
    check_uint64_at(t, 80, 16);
    check_uint64_at(t, 88, 8);
    BOOST_CHECK_EQUAL(t.substr(96, 8), "wibble-2");
    check_uint64_at(t, 104, 32);

    varoom::detail::bytes_input_file_holder i(t);
    varoom::casket c(*i);
    std::vector<std::string> v = c.contents();
    BOOST_REQUIRE_EQUAL(v.size(), 2);
    BOOST_REQUIRE_EQUAL(v[0], "wibble-1");
    BOOST_REQUIRE_EQUAL(v[1], "wibble-2");
    {
        varoom::input_file_holder_ptr ip = c.in("wibble-1");
        std::string u;
        u.resize(16);
        (**ip).read(&u[0], 16);
        BOOST_CHECK_EQUAL(u, s0);
    }
    {
        varoom::input_file_holder_ptr ip = c.in("wibble-2");
        std::string u;
        u.resize(16);
        (**ip).read(&u[0], 16);
        BOOST_CHECK_EQUAL(u, s1);
    }
}
BOOST_AUTO_TEST_CASE( medium_casket_bad )
{
    std::string s0("0123456789abcdef");
    std::string s1("abcdef0123456789");
    std::string t;
    std::function<void(const std::vector<char>&)> f = [&](const std::vector<char>& v) mutable {
        t.clear();
        t.insert(t.end(), v.begin(), v.end());
    };
    {
        varoom::detail::bytes_output_file_holder r(f);
        varoom::casket c(*r);
        {
            varoom::output_file_holder_ptr op = c.out("wibble-1");
            (**op).write(&s0[0], s0.size());
        }
        {
            varoom::output_file_holder_ptr op = c.out("wibble-2");
            (**op).write(&s1[0], s1.size());
        }
    }
    // Expected size:
    //   16     for the string in s0;
    //   16     for the string in s1;
    // + 8      for the toc len;
    // + 8      for the offset of wibble-1
    // + 8      for the length of wibble-1
    // + 8      for the length of the name "wibble-1"
    // + 8      for the characters "wibble-1"
    // + 8      for the offset of wibble-2
    // + 8      for the length of wibble-2
    // + 8      for the length of the name "wibble-2"
    // + 8      for the characters "wibble-2"
    // + 8      for the pointer to the toc
    // = 64
    BOOST_REQUIRE_EQUAL(t.size(), 112);
    BOOST_CHECK_EQUAL(t.substr(0, 16), s0);
    BOOST_CHECK_EQUAL(t.substr(16, 16), s1);
    check_uint64_at(t, 32, 2);
    check_uint64_at(t, 40, 0);
    check_uint64_at(t, 48, 16);
    check_uint64_at(t, 56, 8);
    BOOST_CHECK_EQUAL(t.substr(64, 8), "wibble-1");
    check_uint64_at(t, 72, 16);
    check_uint64_at(t, 80, 16);
    check_uint64_at(t, 88, 8);
    BOOST_CHECK_EQUAL(t.substr(96, 8), "wibble-2");
    check_uint64_at(t, 104, 32);

    varoom::detail::bytes_input_file_holder i(t);
    varoom::casket c(*i);
    std::vector<std::string> v = c.contents();
    BOOST_REQUIRE_EQUAL(v.size(), 2);
    BOOST_REQUIRE_EQUAL(v[0], "wibble-1");
    BOOST_REQUIRE_EQUAL(v[1], "wibble-2");

    varoom::input_file_holder_ptr ip1 = c.in("wibble-1");
    varoom::input_file_holder_ptr ip2 = c.in("wibble-2");
    {
        std::string u;
        u.resize(16);
        (**ip1).read(&u[0], 16);
        BOOST_CHECK_EQUAL(u, s0);
    }
    {
        std::string u;
        u.resize(16);
        (**ip2).read(&u[0], 16);
        BOOST_CHECK_EQUAL(u, s1);
    }
}
