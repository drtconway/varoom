#include "varoom/util/table.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE table tests
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <sstream>

using namespace std;
using namespace varoom;

namespace // anonymous
{
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( streamed_input_table_1 )
{
    using tbl = varoom::detail::streamed_input_table<std::string,uint32_t,uint32_t,std::string>;

    std::string txt =
        "chr1\t11190559\t11190859\tMTOR\n"
        "chr1\t207940883\t207941235\tCD46\n"
        "chr2\t210654202\t210654382\tUNC80\n"
        "chr6\t32165043\t32165403\tNOTCH4\n"
        "chr10\t129908599\t129908839\tMKI67\n"
        "chr11\t118348648\t118348948\tKMT2A\n"
        "chr14\t24709225\t24709405\tTINF2\n"
        "chr16\t89839615\t89839855\tFANCA\n"
        "chr17\t56809784\t56809964\tRAD51C\n"
        "chr21\t44520505\t44520685\tU2AF1\n";

    std::vector<std::string> chrs{"chr1", "chr1", "chr2", "chr6",  "chr10",
                                  "chr11", "chr14", "chr16",  "chr17", "chr21"};
    std::istringstream in(txt);

    tbl t(in);
    tbl::row r;
    while (t.next(r))
    {
        BOOST_REQUIRE_EQUAL(1 <= t.line_number(), true);
        BOOST_REQUIRE_EQUAL(t.line_number() <= chrs.size(), true);
        BOOST_CHECK_EQUAL(std::get<0>(r), chrs[t.line_number() - 1]);
    }
}

BOOST_AUTO_TEST_CASE( streamed_output_table_1 )
{
    using tbl_in  = varoom::detail::streamed_input_table<std::string,uint32_t,uint32_t,std::string>;
    using tbl_out = varoom::detail::streamed_output_table<std::string,uint32_t,uint32_t,std::string>;

    std::string txt =
        "chr1\t11190559\t11190859\tMTOR\n"
        "chr1\t207940883\t207941235\tCD46\n"
        "chr2\t210654202\t210654382\tUNC80\n"
        "chr6\t32165043\t32165403\tNOTCH4\n"
        "chr10\t129908599\t129908839\tMKI67\n"
        "chr11\t118348648\t118348948\tKMT2A\n"
        "chr14\t24709225\t24709405\tTINF2\n"
        "chr16\t89839615\t89839855\tFANCA\n"
        "chr17\t56809784\t56809964\tRAD51C\n"
        "chr21\t44520505\t44520685\tU2AF1\n";

    std::istringstream in(txt);

    tbl_in t(in);
    tbl_in::row r;

    std::ostringstream out;
    tbl_out u(out);

    static_assert(std::is_same<tbl_in::row,tbl_out::row>::value);

    while (t.next(r))
    {
        u << r;
    }

    BOOST_CHECK_EQUAL(out.str(), txt);

}

BOOST_AUTO_TEST_CASE( streamed_output_table_2 )
{
    using tbl_in  = varoom::detail::streamed_input_table<std::string,uint32_t,uint32_t,std::string>;
    using tbl_out = varoom::detail::streamed_output_table<std::string,uint32_t,uint32_t,std::string>;

    std::string txt =
        "chr\tstart\tend\tname\n"
        "chr1\t11190559\t11190859\tMTOR\n"
        "chr1\t207940883\t207941235\tCD46\n"
        "chr2\t210654202\t210654382\tUNC80\n"
        "chr6\t32165043\t32165403\tNOTCH4\n"
        "chr10\t129908599\t129908839\tMKI67\n"
        "chr11\t118348648\t118348948\tKMT2A\n"
        "chr14\t24709225\t24709405\tTINF2\n"
        "chr16\t89839615\t89839855\tFANCA\n"
        "chr17\t56809784\t56809964\tRAD51C\n"
        "chr21\t44520505\t44520685\tU2AF1\n";

    std::istringstream in(txt);

    tbl_in t(in, true);
    tbl_in::row r;

    std::ostringstream out;
    tbl_out u(out, {"chr", "start", "end", "name"});

    static_assert(std::is_same<tbl_in::row,tbl_out::row>::value);

    while (t.next(r))
    {
        u << r;
    }

    BOOST_CHECK_EQUAL(out.str(), txt);

}

BOOST_AUTO_TEST_CASE( inmemory_table_1 )
{
    using tbl_in = varoom::detail::streamed_input_table<std::string,uint32_t,uint32_t,std::string>;
    using tbl_mem = varoom::detail::inmemory_table<std::string,uint32_t,uint32_t,std::string>;

    std::string txt =
        "chr1\t11190559\t11190859\tMTOR\n"
        "chr1\t207940883\t207941235\tCD46\n"
        "chr2\t210654202\t210654382\tUNC80\n"
        "chr6\t32165043\t32165403\tNOTCH4\n"
        "chr10\t129908599\t129908839\tMKI67\n"
        "chr11\t118348648\t118348948\tKMT2A\n"
        "chr14\t24709225\t24709405\tTINF2\n"
        "chr16\t89839615\t89839855\tFANCA\n"
        "chr17\t56809784\t56809964\tRAD51C\n"
        "chr21\t44520505\t44520685\tU2AF1\n";

    std::istringstream in(txt);

    tbl_in t(in);
    tbl_in::row r;

    tbl_mem u;

    while (t.next(r))
    {
        u.push_back(r);
    }
}

BOOST_AUTO_TEST_CASE( table_map_1 )
{
    using tbl_in = varoom::detail::streamed_input_table<std::string,uint32_t,uint32_t,std::string>;
    using tbl_mem = varoom::detail::inmemory_table<std::string,uint32_t,uint32_t,std::string,uint32_t>;
    using tbl_out = varoom::detail::streamed_outmemory_table<std::string,uint32_t,uint32_t,std::string,uint32_t>;

    std::string txt =
        "chr1\t11190559\t11190859\tMTOR\n"
        "chr1\t207940883\t207941235\tCD46\n"
        "chr2\t210654202\t210654382\tUNC80\n"
        "chr6\t32165043\t32165403\tNOTCH4\n"
        "chr10\t129908599\t129908839\tMKI67\n"
        "chr11\t118348648\t118348948\tKMT2A\n"
        "chr14\t24709225\t24709405\tTINF2\n"
        "chr16\t89839615\t89839855\tFANCA\n"
        "chr17\t56809784\t56809964\tRAD51C\n"
        "chr21\t44520505\t44520685\tU2AF1\n";

    std::istringstream in(txt);

    tbl_in t(in);

    tbl_mem u;
    tbl_out v(u);

    std::function<void(const tbl_in::row&, tbl_out::row&)> f = [](const tbl_in::row& p_x, tbl_out::row& p_y) {
        varoom::table_utils::copy<0,4>(p_x, p_y);
        std::get<4>(p_y) = std::get<2>(p_x) - std::get<1>(p_x);
    };

    varoom::table_utils::map(f, t, v);

    BOOST_REQUIRE_EQUAL(u.size(), 10);
}
