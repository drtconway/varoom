#ifndef VAROOM_UTIL_GAMMA_ESTIMATOR_HPP
#define VAROOM_UTIL_GAMMA_ESTIMATOR_HPP

#include <cmath>
#include <utility>

namespace varoom
{
    class gamma_estimator
    {
    public:
        gamma_estimator()
            : m_n(0), m_sx(0), m_slx(0), m_sxlx(0), m_last(1e-6)
        {
        }

        gamma_estimator& push_back(const double& p_x)
        {
            ++m_n;
            if (p_x > 0)
            {
                double lx = std::log(p_x);
                m_sx += p_x;
                m_slx += lx;
                m_sxlx += p_x*lx;
                m_last = p_x;
            }
            return *this;
        }

        size_t count() const
        {
            return m_n;
        }

        std::pair<double,double> operator()() const
        {
            size_t n = m_n;
            double sx = m_sx;
            double slx = m_slx;
            double sxlx = m_sxlx;
            double v = n * sxlx - sx*slx;
            if (v == 0.0)
            {
                // If all the observations were equal,
                // v comes out as 0. In that case add
                // a small jitter observation to avoid
                // the divide-by-zero.
                double x = 1.000001 * m_last;
                double lx = std::log(x);
                n += 1;
                sx += x;
                slx += lx;
                sxlx += x*lx;
                v = n * sxlx - sx*slx;
            }
            double k = n * sx / v;
            double t = v / (n*n);
            return std::pair<double,double>(k, t);
        }

    private:
        size_t m_n;
        double m_sx;
        double m_slx;
        double m_sxlx;
        double m_last;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_GAMMA_ESTIMATOR_HPP
