#ifndef VAROOM_UTIL_GAMMA_ESTIMATOR_HPP
#define VAROOM_UTIL_GAMMA_ESTIMATOR_HPP

#include <cmath>

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

        std::pair<double,double> operator()() const
        {
            double v = m_n * m_sxlx - m_sx*m_slx;
            double k = m_n * m_sx / v;
            double t = v / (m_n*m_n);
            return std::pair<double,double>(k, t);
        }

    private:
        size_t m_n;
        double m_sx;
        double m_slx;
        double m_sxlx;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_GAMMA_ESTIMATOR_HPP
