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
            : m_n(0), m_sx(0), m_slx(0), m_sxlx(0)
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
            if (v > 1e-12)
            {
                double k = n * sx / v;
                double t = v / (n*n);
                return std::pair<double,double>(k, t);
            }
            else
            {
                // Make the variance 1% of the mean.
                //
                // use mean = k*th, variance = k*th*th.
                //
                // variance = mean * 0.01 => th = 0.01
                double m = sx / double(n);
                double t = 0.01;
                double k = m / t;
                return std::pair<double,double>(k, t);
            }
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
