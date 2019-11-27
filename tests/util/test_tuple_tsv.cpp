#include "varoom/util/tuple_tsv.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tuple_tsv tests
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <tuple>
#include <boost/any.hpp>

namespace // anonymous
{
    std::string conv_str(const varoom::subtext& p_s)
    {
        return p_s;
    }
    uint32_t conv_uint32_t(const varoom::subtext& p_s)
    {
        auto r = boost::make_iterator_range(p_s.first, p_s.second);
        return boost::lexical_cast<uint32_t>(r);
    }
    size_t conv_size_t(const varoom::subtext& p_s)
    {
        auto r = boost::make_iterator_range(p_s.first, p_s.second);
        return boost::lexical_cast<size_t>(r);
    }

    std::pair<std::string,uint32_t> conv_locus(const varoom::subtext& p_s)
    {
        std::vector<varoom::subtext> parts;
        p_s.split(':', parts);
        auto r = boost::make_iterator_range(parts[1].first, parts[1].second);
        std::pair<std::string,uint32_t> v(parts[0], boost::lexical_cast<uint32_t>(r));
        return v;
    }
}
// namespace anonymous

using namespace std;
using namespace varoom;

BOOST_AUTO_TEST_CASE( cast1 )
{
    string txt = "chr1\t123\t234\t345\t456\t567";
    vector<subtext> ss;
    subtext(txt).split('\t', ss);
    tuple<string,uint32_t,size_t,size_t,size_t,size_t> t;
    tuple_tsv::to_tuple(ss, t);
    BOOST_CHECK_EQUAL(std::get<0>(t), "chr1");
    BOOST_CHECK_EQUAL(std::get<1>(t), 123);
    BOOST_CHECK_EQUAL(std::get<2>(t), 234);
    BOOST_CHECK_EQUAL(std::get<3>(t), 345);
    BOOST_CHECK_EQUAL(std::get<4>(t), 456);
    BOOST_CHECK_EQUAL(std::get<5>(t), 567);
}

BOOST_AUTO_TEST_CASE( cast2 )
{
    std::function<std::string(const subtext&)> f = conv_str;
    std::function<uint32_t(const subtext&)> g = conv_uint32_t;
    std::function<size_t(const subtext&)> h = conv_size_t;

    string txt = "chr1\t123\t234\t345\t456\t567";
    vector<boost::any> ts{tuple_tsv::to_str(), tuple_tsv::to_uint32(),
                          tuple_tsv::to_size_t(), tuple_tsv::to_size_t(),
                          tuple_tsv::to_size_t(), tuple_tsv::to_size_t()};
    vector<subtext> ss;
    subtext(txt).split('\t', ss);
    tuple<string,uint32_t,size_t,size_t,size_t,size_t> t;
    tuple_tsv::to_tuple(ss, t, ts);
    BOOST_CHECK_EQUAL(std::get<0>(t), "chr1");
    BOOST_CHECK_EQUAL(std::get<1>(t), 123);
    BOOST_CHECK_EQUAL(std::get<2>(t), 234);
    BOOST_CHECK_EQUAL(std::get<3>(t), 345);
    BOOST_CHECK_EQUAL(std::get<4>(t), 456);
    BOOST_CHECK_EQUAL(std::get<5>(t), 567);
}

BOOST_AUTO_TEST_CASE( cast3 )
{
    std::function<std::pair<std::string,uint32_t>(const subtext&)> f = conv_locus;
    std::function<size_t(const subtext&)> h = conv_size_t;

    string txt = "chr1:123\t234\t345\t456\t567";
    vector<boost::any> ts{f, h, h, h, h};
    vector<subtext> ss;
    subtext(txt).split('\t', ss);
    tuple<std::pair<string,uint32_t>,size_t,size_t,size_t,size_t> t;
    tuple_tsv::to_tuple(ss, t, ts);
    BOOST_CHECK_EQUAL(std::get<0>(t).first, "chr1");
    BOOST_CHECK_EQUAL(std::get<0>(t).second, 123);
    BOOST_CHECK_EQUAL(std::get<1>(t), 234);
    BOOST_CHECK_EQUAL(std::get<2>(t), 345);
    BOOST_CHECK_EQUAL(std::get<3>(t), 456);
    BOOST_CHECK_EQUAL(std::get<4>(t), 567);
}
