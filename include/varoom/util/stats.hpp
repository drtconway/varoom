#ifndef VAROOM_UTIL_STATS_HPP
#define VAROOM_UTIL_STATS_HPP

#include <math>
#include <utility>
#include <stdexcept>
#include <vector>

namespace varoom
{
    namespace stats
    {
        double kl_divergence(const std::vector<double>& p_ps, const std::vector<double>& p_qs)
        {
            if (p_ps.size() != p_qs.size())
            {
                throw std::invalid_argument("both distributions must have the same number of categories");
            }

            double d = 0;
            for (size_t i = 0; i < p_ps.size(); ++i)
            {
                double p = p_ps[i];
                double q = p_qs[i];
                if (p == 0.0)
                {
                    continue;
                }
                if (q == 0.0)
                {
                    throw std::domain_error("Q probabilities must be strictly > 0");
                }
                d += p*std::log(q/p);
            }
            return d;
        }

        std::pair<double,double> gamma_mm(const std::vector<double>& p_xs)
        {
            double sx = 0;
            double sxl = 0;
            double sl = 0;
            double n = p_xs.size();
            for (size_t i = 0; i < p_xs.size(); ++i)
            {
                double x = p_xs[i];
                double lx = std::log(x);
                sx += x;
                sxl += x*lx;
                sl += lx;
            }
            double k = n*sx / (n*sxl - sl*sx);
            double theta = (n*sxl - sl*sx) / (n*n);
            return std::pair<double,double>(k, theta);
        }
    }
    // namespace stats
}
// namespace varoom

#endif // VAROOM_UTIL_STATS_HPP
