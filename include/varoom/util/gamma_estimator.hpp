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

        gamma_estimator(const size_t& p_n, const double& p_sx, const double& p_slx, const double& p_sxlx)
            : m_n(p_n), m_sx(p_sx), m_slx(p_slx), m_sxlx(p_sxlx)
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

        double sum_x() const
        {
            return m_sx;
        }

        double sum_log_x() const
        {
            return m_slx;
        }

        double sum_x_log_x() const
        {
            return m_sxlx;
        }

        gamma_estimator& operator+=(const gamma_estimator& p_other)
        {
            m_n += p_other.m_n;
            m_sx += p_other.m_sx;
            m_slx += p_other.m_slx;
            m_sxlx += p_other.m_sxlx;
            return *this;
        }

        std::pair<double,double> operator()() const
        {
            return estimate(m_n, m_sx, m_slx, m_sxlx);
        }

        static std::pair<double,double> estimate(const size_t& p_n, const double& p_sx, const double& p_slx, const double& p_sxlx)
        {
            size_t n = p_n;
            double sx = p_sx;
            double slx = p_slx;
            double sxlx = p_sxlx;
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
