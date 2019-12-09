#include "varoom/util/concurrent_deque.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE concurrent_deque tests
#include <boost/test/unit_test.hpp>

#include <set>
#include <thread>
#include <iostream>

namespace // anonymous
{
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( enc1 )
{
    const int N = 1000;
    varoom::concurrent_deque<int> q(10);
    std::thread c([&]() {
        int x = 0;
        std::set<int> wanted;
        for (int i = 0; i < N; ++i)
        {
            wanted.insert(i);
        }
        while (q.pop_front(x))
        {
            BOOST_CHECK_EQUAL(wanted.count(x), true);
            wanted.erase(x);
        }
        BOOST_CHECK_EQUAL(wanted.size(), 0);
    });
    std::thread p([&]() {
        for (int i = 0; i < N; ++i)
        {
            q.push_back(i);
        }
        q.end();
    });

    p.join();
    c.join();
    BOOST_CHECK_EQUAL(q.size(), 0);
}
