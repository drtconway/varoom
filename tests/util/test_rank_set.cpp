#include "varoom/util/rank_set.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE rank_set tests
#include <boost/test/unit_test.hpp>

#include <random>
#include <set>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/range/iterator_range.hpp>

namespace // anonymous
{
    using tree_ptr = varoom::detail::avl_tree_ptr<int>;
    using stree = varoom::detail::avl_tree<std::string>;
    using stree_ptr = varoom::detail::avl_tree_ptr<std::string>;

    void traverse(tree_ptr x, std::string& r)
    {
        if (!x)
        {
            r += ".";
            return;
        }
        r += "(";
        traverse(x->m_lhs, r);
        r += boost::lexical_cast<std::string>(x->m_key);
        traverse(x->m_rhs, r);
        r += ")";
    }
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( rot_1 )
{
    std::string v;
    std::function<void(const int&)> f = [&] (const int& p_x) mutable {
        v.push_back(p_x);
    };

    stree_ptr x(new stree("x"));
    stree_ptr y(new stree("y"));
    stree_ptr t1(new stree("t1"));
    stree_ptr t2(new stree("t2"));
    stree_ptr t3(new stree("t3"));
    
    x->m_lhs = t1;
    x->m_rhs = t2;
    x->refresh_size_and_height();
    BOOST_CHECK_EQUAL(x->size(), 3);
    BOOST_CHECK_EQUAL(x->height(), 2);
    BOOST_CHECK_EQUAL(x->balance_factor(), 0);

    y->m_lhs = x;
    y->m_rhs = t3;
    y->refresh_size_and_height();

    BOOST_CHECK_EQUAL(y->size(), 5);
    BOOST_CHECK_EQUAL(y->height(), 3);
    BOOST_CHECK_EQUAL(y->balance_factor(), 1);

    stree_ptr z;
    z = varoom::detail::avl_tree<std::string>::rotate_right(y);
    BOOST_CHECK_EQUAL(z.get(), x.get());
    BOOST_CHECK_EQUAL(x->size(), 5);
    BOOST_CHECK_EQUAL(x->height(), 3);
    BOOST_CHECK_EQUAL(x->balance_factor(), -1);
    BOOST_CHECK_EQUAL(y->size(), 3);
    BOOST_CHECK_EQUAL(y->height(), 2);
    BOOST_CHECK_EQUAL(y->balance_factor(), 0);

    z = stree::rotate_left(x);
    BOOST_CHECK_EQUAL(z.get(), y.get());
    BOOST_CHECK_EQUAL(x->size(), 3);
    BOOST_CHECK_EQUAL(x->height(), 2);
    BOOST_CHECK_EQUAL(x->balance_factor(), 0);
    BOOST_CHECK_EQUAL(y->size(), 5);
    BOOST_CHECK_EQUAL(y->height(), 3);
    BOOST_CHECK_EQUAL(y->balance_factor(), 1);
}

BOOST_AUTO_TEST_CASE( avl_1 )
{
    std::vector<int> v;
    std::function<void(const int&)> f = [&] (const int& p_x) mutable {
        v.push_back(p_x);
    };

    varoom::detail::avl_tree_ptr<int> x;
    {
        v.clear();
        varoom::detail::avl_tree<int>::visit(x, f);
        BOOST_CHECK_EQUAL(v.size(), 0);
    }
    x = varoom::detail::avl_tree<int>::insert(1, x);
    {
        v.clear();
        varoom::detail::avl_tree<int>::visit(x, f);
        BOOST_REQUIRE_EQUAL(v.size(), 1);
        BOOST_REQUIRE_EQUAL(v[0], 1);
    }
    x = varoom::detail::avl_tree<int>::insert(2, x);
    {
        v.clear();
        varoom::detail::avl_tree<int>::visit(x, f);
        BOOST_REQUIRE_EQUAL(v.size(), 2);
        BOOST_REQUIRE_EQUAL(v[0], 1);
        BOOST_REQUIRE_EQUAL(v[1], 2);
    }
    x = varoom::detail::avl_tree<int>::insert(3, x);
    {
        v.clear();
        varoom::detail::avl_tree<int>::visit(x, f);
        BOOST_REQUIRE_EQUAL(v.size(), 3);
        BOOST_REQUIRE_EQUAL(v[0], 1);
        BOOST_REQUIRE_EQUAL(v[1], 2);
        BOOST_REQUIRE_EQUAL(v[2], 3);
    }

    std::string r;
    traverse(x, r);
    BOOST_CHECK_EQUAL(r, "((.1.)2(.3.))");
}

BOOST_AUTO_TEST_CASE( avl_2 )
{
    std::vector<int> v;
    std::function<void(const int&)> f = [&] (const int& p_x) mutable {
        v.push_back(p_x);
    };

    varoom::detail::avl_tree_ptr<int> x;
    {
        v.clear();
        varoom::detail::avl_tree<int>::visit(x, f);
        BOOST_CHECK_EQUAL(v.size(), 0);
    }
    x = varoom::detail::avl_tree<int>::insert(2, x);
    {
        v.clear();
        varoom::detail::avl_tree<int>::visit(x, f);
        BOOST_REQUIRE_EQUAL(v.size(), 1);
        BOOST_REQUIRE_EQUAL(v[0], 2);
    }
    x = varoom::detail::avl_tree<int>::insert(0, x);
    {
        v.clear();
        varoom::detail::avl_tree<int>::visit(x, f);
        BOOST_REQUIRE_EQUAL(v.size(), 2);
        BOOST_REQUIRE_EQUAL(v[0], 0);
        BOOST_REQUIRE_EQUAL(v[1], 2);
    }
    x = varoom::detail::avl_tree<int>::insert(3, x);
    {
        v.clear();
        varoom::detail::avl_tree<int>::visit(x, f);
        BOOST_REQUIRE_EQUAL(v.size(), 3);
        BOOST_REQUIRE_EQUAL(v[0], 0);
        BOOST_REQUIRE_EQUAL(v[1], 2);
        BOOST_REQUIRE_EQUAL(v[2], 3);
    }
    x = varoom::detail::avl_tree<int>::insert(1, x);
    {
        v.clear();
        varoom::detail::avl_tree<int>::visit(x, f);
        BOOST_REQUIRE_EQUAL(v.size(), 4);
        BOOST_REQUIRE_EQUAL(v[0], 0);
        BOOST_REQUIRE_EQUAL(v[1], 1);
        BOOST_REQUIRE_EQUAL(v[2], 2);
        BOOST_REQUIRE_EQUAL(v[3], 3);
    }
}

BOOST_AUTO_TEST_CASE( avl_3 )
{
    std::vector<int> v;
    std::function<void(const int&)> f = [&] (const int& p_x) mutable {
        v.push_back(p_x);
    };

    varoom::detail::avl_tree_ptr<int> x;
    size_t N = 1000;
    std::set<int> S;
    std::mt19937_64 rng(17);
    std::uniform_int_distribution<uint64_t> U(0, 1<<30);
    for (size_t i = 0; i < N; ++i)
    {
        int a = U(rng);
        //std::cerr << "inserting: " << a << std::endl;
        S.insert(a);
        x = varoom::detail::avl_tree<int>::insert(a, x);

        v.clear();
        varoom::detail::avl_tree<int>::visit(x, f);
        BOOST_REQUIRE_EQUAL(v.size(), S.size());
        for (size_t j = 1; j < v.size(); ++j)
        {
            BOOST_REQUIRE_EQUAL(v[j-1] < v[j], true);
        }

        std::string r;
        traverse(x, r);
        //std::cerr << i << '\t' << r << std::endl;
    }
    BOOST_CHECK_EQUAL(x->size(), S.size());

    v.clear();
    varoom::detail::avl_tree<int>::visit(x, f);
    BOOST_REQUIRE_EQUAL(v.size(), S.size());
    size_t i = 0;
    for (auto s = S.begin(); s != S.end(); ++s, ++i)
    {
        BOOST_CHECK_EQUAL(x->contains(*s), true);
        BOOST_CHECK_EQUAL(x->rank(*s), i);
        BOOST_CHECK_EQUAL(x->select(i), *s);
        BOOST_CHECK_EQUAL(v[i], *s);
    }

    std::shuffle(v.begin(), v.end(), rng);
    std::vector<int> w(v.begin(), v.end());
    for (size_t i = 0; i < v.size(); ++i)
    {
        int a = v[i];
        BOOST_CHECK_EQUAL(S.count(a), true);
        BOOST_CHECK_EQUAL(x->contains(a), true);
        //std::cerr << i << '\t' << a << std::endl;
        x = varoom::detail::avl_tree<int>::remove(a, x);
        varoom::detail::avl_tree<int>::sane(x);
        S.erase(a);
        BOOST_CHECK_EQUAL((x ? x->size() : 0), S.size());
        BOOST_CHECK_EQUAL(S.count(a), false);
        if (x)
        {
            BOOST_CHECK_EQUAL(x->contains(a), false);
        }
    }
}

BOOST_AUTO_TEST_CASE( rank_set_1 )
{
    varoom::rank_set<size_t> R;
    R.insert(6);
    R.insert(1);
    BOOST_CHECK_EQUAL(R.rank(8), 2);
}
