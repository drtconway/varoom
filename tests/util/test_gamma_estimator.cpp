#include "varoom/util/gamma_estimator.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE gamma_estimator tests
#include <boost/test/unit_test.hpp>

#include <iostream>

using namespace std;
using namespace varoom;

namespace // anonymous
{
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( enc1 )
{
    vector<double> xs{0.0240453, 0.0240453};

    gamma_estimator g;
    for (size_t i = 0; i < xs.size(); ++i)
    {
        g.push_back(xs[i]);
    }

    pair<double,double> kth = g();

    BOOST_CHECK_EQUAL(std::isfinite(kth.first), true);
    BOOST_CHECK_EQUAL(std::isfinite(kth.second), true);
}

BOOST_AUTO_TEST_CASE( enc2 )
{
    vector<double> xs;
    size_t N = 100;
    for (size_t i = 0; i < N; ++i)
    {
        xs.push_back(0.0240453);
    }

    gamma_estimator g;
    for (size_t i = 0; i < xs.size(); ++i)
    {
        g.push_back(xs[i]);
    }

    pair<double,double> kth = g();

    //cout << kth.first << '\t' << kth.second << endl;
    BOOST_CHECK_EQUAL(std::isfinite(kth.first), true);
    BOOST_CHECK_EQUAL(std::isfinite(kth.second), true);
}
