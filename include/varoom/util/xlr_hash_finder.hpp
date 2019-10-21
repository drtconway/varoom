#ifndef VAROOM_XLR_HASH_FINDER_HPP
#define VAROOM_XLR_HASH_FINDER_HPP

#ifndef VAROOM_XLR_HASH_HPP
#include "varoom/util/xlr_hash.hpp"
#endif

#include <random>
#include <set>
#include <sstream>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/miller_rabin.hpp>

#define BOOST_LOG_DYN_LINK 1
#include <boost/log/trivial.hpp>

namespace varoom
{
    namespace detail
    {
        class bit_counter
        {
        public:
            bit_counter()
                : m_counts(64), m_n(0)
            {
            }

            void add(std::uint64_t p_x)
            {
                for (size_t i = 0; i < 64; ++i)
                {
                    m_counts[i] += p_x & 1;
                    p_x >>= 1;
                }
                m_n += 1;
            }

            double chi2() const
            {
                double m = 0.5 * m_n;

                double v = 0;
                for (size_t i = 0; i < 64; ++i)
                {
                    double w = m_counts[i] - m;
                    v += w*w / m;
                }
                return v;
            }

        private:
            std::vector<std::uint64_t> m_counts;
            std::uint64_t m_n;
        };

        template <int B>
        double chi2(const std::uint64_t& p_x)
        {
            const std::uint64_t N = 1 << B;
            const std::uint64_t M = (N - 1);
            const std::uint64_t J = 64 / B;
            size_t cs[N];
            for (size_t i = 0; i < N; ++i)
            {
                cs[i] = 0;
            }
            std::uint64_t x = p_x;
            for (size_t i = 0; i < J; ++i)
            {
                cs[x&M] += 1;
                x >>= B;
            }

            const double m = double(J) / double(N);
            double c2 = 0.0;
            for (size_t i = 0; i < N; ++i)
            {
                double d = cs[i] - m;
                c2 += d*d / m;
            }
            return c2;
        }

        std::string as_binary(const std::uint64_t& p_x)
        {
            std::string s;
            std::uint64_t x = p_x;
            for (size_t i = 0; i < 64; ++i)
            {
                s.push_back("01"[x&1]);
                x >>= 1;
            }
            return std::string(s.rbegin(), s.rend());
        }

        std::string as_hex(const std::uint64_t& p_x)
        {
            std::string s;
            std::uint64_t x = p_x;
            for (size_t i = 0; i < 16; ++i)
            {
                s.push_back("0123456789abcdef"[x&15]);
                x >>= 4;
            }
            return std::string(s.rbegin(), s.rend());
        }
        std::string as_json(const std::vector<std::uint64_t>& p_xs)
        {
            std::ostringstream out;
            out << "{";
            for (size_t i = 0; i < p_xs.size(); ++i)
            {
                if (i > 0)
                {
                    out << ", ";
                }
                out << p_xs[i];
            }
            out << "}";
            return out.str();
        }
    }
    // namespace detail

    class xlr_hash_finder
    {
    public:
        xlr_hash_finder(const std::uint64_t& p_num_blocks,
                        const std::uint64_t& p_seed,
                        const std::vector<std::string>& p_strings)
            : m_num_blocks(p_num_blocks), m_strings(p_strings), m_rng(p_seed), m_primes64(65536, 0)
        {
        }

        void find(const std::uint64_t& p_num_trials, std::vector<std::uint64_t>& p_res)
        {
            std::vector<std::uint64_t> best;
            std::vector<double> bestScores;
            make_coeffs(best);
            evaluate(best, bestScores);

            std::vector<std::uint64_t> next;
            std::vector<double> nextScores;
            for (size_t i = 1; i < p_num_trials; ++i)
            {
                next.clear();
                nextScores.clear();

                make_coeffs(next);
                evaluate(next, nextScores);

                if (is_better(nextScores, bestScores))
                {
                    next.swap(best);
                    nextScores.swap(bestScores);
                }
                if ((i & 0xf) == 0)
                {
                    BOOST_LOG_TRIVIAL(info) << i << '\t' << detail::as_json(best);
                    for (size_t j = 0; j < bestScores.size(); ++j)
                    {
                        BOOST_LOG_TRIVIAL(info) << i << '\t' << j << '\t' << bestScores[j];
                    }
                }
            }
            BOOST_LOG_TRIVIAL(info) << "final: " << detail::as_json(best);
            for (size_t j = 0; j < bestScores.size(); ++j)
            {
                BOOST_LOG_TRIVIAL(info) << "final: " << j << '\t' << bestScores[j];
            }
            p_res.swap(best);
        }

    private:
        void make_coeffs(std::vector<std::uint64_t>& p_coeffs)
        {
            p_coeffs.clear();

            std::uniform_int_distribution<std::uint64_t> X(0, m_primes64.size() - 1);
            std::uniform_int_distribution<std::uint64_t> U(7, 31);
            for (size_t i = 0; i < m_num_blocks; ++i)
            {
                p_coeffs.push_back(getPrime(X(m_rng)));
                p_coeffs.push_back(U(m_rng));
                p_coeffs.push_back(U(m_rng));
            }
        }

        void evaluate(const std::vector<std::uint64_t>& p_coeffs,
                      std::vector<double>& p_scores)
        {
            xlr_hash H(p_coeffs);
            p_scores.clear();

            const std::uint64_t S = 0x317371dc620c3879ULL;

            // Test the first million integers with a zero seed.
            //
            {
                detail::bit_counter b;
                for (std::uint64_t x = 0; x < 1024*1024; ++x)
                {
                    std::uint64_t h = H(0, x);
                    b.add(h);
                }
                p_scores.push_back(b.chi2());
            }

            // Test the first million integers with a big prime seed.
            //
            {
                detail::bit_counter b;
                for (std::uint64_t x = 0; x < 1024*1024; ++x)
                {
                    std::uint64_t h = H(S, x);
                    b.add(h);
                }
                p_scores.push_back(b.chi2());
            }

            // Test a million random integers with a zero seed.
            //
            {
                std::uniform_int_distribution<std::uint64_t> U(0, 1ULL << 63);
                detail::bit_counter b;
                for (size_t i = 0; i < 1024*1024; ++i)
                {
                    std::uint64_t x = U(m_rng);
                    std::uint64_t h = H(0, x);
                    b.add(h);
                }
                p_scores.push_back(b.chi2());
            }

            // Test a million random integers with a big prime seed.
            //
            {
                std::uniform_int_distribution<std::uint64_t> U(0, 1ULL << 63);
                detail::bit_counter b;
                for (size_t i = 0; i < 1024*1024; ++i)
                {
                    std::uint64_t x = U(m_rng);
                    std::uint64_t h = H(S, x);
                    b.add(h);
                }
                p_scores.push_back(b.chi2());
            }

            // Test the 1-bit avalanche for 1024 random integers.
            //
            {
                std::uniform_int_distribution<std::uint64_t> U(0, 1ULL << 63);
                detail::bit_counter b;
                for (size_t i = 0; i < 1024; ++i)
                {
                    std::uint64_t x = U(m_rng);
                    std::uint64_t h0 = H(S, x);
                    for (size_t j = 0; j < 64; ++j)
                    {
                        std::uint64_t y = x ^ (1ULL << j);
                        std::uint64_t h1 = H(S, y);
                        if (0)
                        {
                            BOOST_LOG_TRIVIAL(debug) << detail::as_hex(x)
                                                     << '\t' << detail::as_hex(y)
                                                     << '\t' << detail::as_hex(x ^ y)
                                                     << '\t' << detail::as_binary(h0 ^ h1)
                                                     ;
                        }
                        b.add(h0 ^ h1);
                    }
                }
                p_scores.push_back(b.chi2());
            }

            // Test the strings with a zero seed.
            //
            {
                detail::bit_counter b;
                for (size_t i = 0; i < m_strings.size(); ++i)
                {
                    std::uint64_t h = H(0, m_strings[i]);
                    b.add(h);
                }
                p_scores.push_back(b.chi2());
            }

            // Test the strings with a big prime seed.
            //
            {
                detail::bit_counter b;
                for (size_t i = 0; i < m_strings.size(); ++i)
                {
                    std::uint64_t h = H(S, m_strings[i]);
                    b.add(h);
                }
                p_scores.push_back(b.chi2());
            }

            if (0)
            {
                for (size_t i = 0; i < p_scores.size(); ++i)
                {
                    BOOST_LOG_TRIVIAL(debug) << i << '\t' << p_scores[i];
                }
            }
        }

        static bool is_better(const std::vector<double>& p_lhs, const std::vector<double>& p_rhs)
        {
            double l = 0;
            double r = 0;
            for (size_t i = 0; i < p_lhs.size(); ++i)
            {
                l += p_lhs[i];
                r += p_rhs[i];
            }
            return l < r;
        }

        std::uint64_t getPrime(size_t p_idx)
        {
            if (m_primes64[p_idx] == 0)
            {
                m_primes64[p_idx] = bigPrime(64);
            }
            return m_primes64[p_idx];
        }

        std::uint64_t bigPrime(size_t p_num_bits)
        {
            std::uniform_int_distribution<std::uint64_t> U(0, (1ULL << p_num_bits) - 1);

            while (true)
            {
                std::uint64_t u = U(m_rng) | 1ULL;
                boost::multiprecision::cpp_int n = u;
                if (miller_rabin_test(n, 25, m_rng) && miller_rabin_test((n-1)/2, 25, m_rng))
                {
                    return u;
                }
            }
        }

        const std::uint64_t m_num_blocks;
        const std::vector<std::string>& m_strings;
        std::mt19937_64 m_rng;

        std::vector<std::uint64_t> m_primes64;
    };
}
// namespace

#endif // VAROOM_XLR_HASH_FINDER_HPP
