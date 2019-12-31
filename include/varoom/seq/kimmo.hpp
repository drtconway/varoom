#ifndef VAROOM_SEQ_KIMMO_HPP
#define VAROOM_SEQ_KIMMO_HPP

#include <algorithm>
#include <array>
#include <string>
#include <unordered_map>
#include <vector>
#include <nlohmann/json.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>

#ifndef VAROOM_UTIL_FILES_HPP
#include "varoom/util/files.hpp"
#endif

#include <chrono>
#include <cmath>
#include <iostream>

namespace varoom
{
    namespace seq
    {
        class kimmo
        {
        public:
            kimmo(const std::string& p_basename)
            {
                {
                    varoom::input_file_holder_ptr inp = varoom::files::in(p_basename + ".toc");
                    nlohmann::json jj;
                    (**inp) >> jj;
                    m_size = jj["size"];
                    {
                        nlohmann::json kk = jj["codes"];
                        for (size_t i = 0; i < kk.size(); ++i)
                        {
                            std::vector<char> v = kk[i];
                            m_codes.push_back(v);
                        }
                    }
                    {
                        nlohmann::json kk = jj["index"];
                        for (size_t i = 0; i < kk.size(); ++i)
                        {
                            uint8_t c = kk[i][0];
                            size_t n = kk[i][1];
                            size_t j = kk[i][2];
                            m_code_index[c] = std::pair<size_t,size_t>(n, j);
                        }
                    }
                }
                sdsl::load_from_file(m_S, p_basename + "-S.bv");
                sdsl::load_from_file(m_D, p_basename + "-D.bv");
                m_select_1_D = sdsl::bit_vector::select_1_type(&m_D);
            }

            size_t size() const
            {
                return m_size;
            }

            char operator[](size_t p_idx) const
            {
                size_t p0 = m_select_1_D(p_idx+1);
                size_t p1 = m_select_1_D(p_idx+2);
                size_t n = p1 - p0;
                size_t j = 0;
                for (size_t p = p0; p < p1; ++p)
                {
                    j = (j << 1) | m_S[p];
                }
                return m_codes[n-1][j];
            }

            std::string slice(size_t p_begin, size_t p_end) const
            {
                size_t pb = m_select_1_D(p_begin+1);
                size_t pe = m_select_1_D(p_end+1);
                std::string s;
                s.reserve(p_end - p_begin);
                size_t n = 0;
                size_t j = 0;
                for (size_t p = pb; p < pe; ++p)
                {
                    if (m_D[p] == 1)
                    {
                        if (n > 0)
                        {
                            s.push_back(m_codes[n-1][j]);
                            n = 0;
                            j = 0;
                        }
                    }
                    n += 1;
                    j = (j << 1) | m_S[p];
                }
                if (n > 0)
                {
                    s.push_back(m_codes[n-1][j]);
                }
                return s;
            }

            static void make(const std::string p_text, const std::string& p_basename)
            {
                std::array<size_t,256> cts;
                for (size_t i = 0; i < cts.size(); ++i)
                {
                    cts[i] = 0;
                }

                for (size_t i = 0; i < p_text.size(); ++i)
                {
                    cts[p_text[i]] += 1;
                }

                using count_and_char = std::pair<size_t,uint8_t>;

                std::vector<count_and_char> items;
                for (size_t i = 0; i < cts.size(); ++i)
                {
                    if (cts[i] == 0)
                    {
                        continue;
                    }
                    items.push_back(count_and_char(cts[i], uint8_t(i)));
                }
                std::sort(items.rbegin(), items.rend());

                std::vector<std::vector<char>> codes;
                std::unordered_map<char,std::pair<size_t,size_t>> code_index;

                size_t n = 1;
                size_t j = 0;
                size_t total_bits = 0;
                for (size_t i = 0; i < items.size(); ++i)
                {
                    if (codes.size() < n)
                    {
                        codes.push_back(std::vector<char>(1ULL << n, '\0'));
                    }
                    codes[n-1][j] = items[i].second;
                    code_index[items[i].second] = std::pair<size_t,size_t>(n, j);

                    total_bits += n * items[i].first;

                    std::cout << i << '\t' << n << '\t' << j << '\t' << char(items[i].second) << '\t' << items[i].first << std::endl;

                    j += 1;
                    if (j == (1ULL << n))
                    {
                        n += 1;
                        j = 0;
                    }
                }
                //std::cout << "total_bits = " << total_bits << std::endl;
                {
                    double e = 0.0;
                    for (size_t i = 0; i < items.size(); ++i)
                    {
                        double p = double(items[i].first)/double(p_text.size());
                        e += -p*std::log(p);
                    }
                    e /= std::log(2);
                    std::cerr << "E0 = " << e << std::endl;
                }

                {
                    varoom::output_file_holder_ptr outp = varoom::files::out(p_basename + ".toc");
                    nlohmann::json jj = nlohmann::json::object();
                    jj["size"] = p_text.size();
                    jj["codes"] = codes;
                    jj["index"] = nlohmann::json::array();
                    for (auto itr = code_index.begin(); itr != code_index.end(); ++itr)
                    {
                        jj["index"].push_back({itr->first, itr->second.first, itr->second.second});
                    }
                    (**outp) << jj << std::endl;
                }

                sdsl::bit_vector S(total_bits, 0);
                sdsl::bit_vector D(total_bits+1, 0);

                size_t b = 0;
                for (size_t i = 0; i < p_text.size(); ++i)
                {
                    auto itr = code_index.find(p_text[i]);
                    size_t n = itr->second.first;
                    size_t j = itr->second.second;
                    for (size_t k = 0; k < n; ++k, ++b)
                    {
                        size_t w = (n - 1) - k;
                        S[b] = (j >> w) & 1;
                        D[b] = (k == 0);
                    }
                }
                D[b] = 1;

                sdsl::store_to_file(S, p_basename + "-S.bv");
                sdsl::store_to_file(D, p_basename + "-D.bv");

                double p1 = double(p_text.size()) / double(total_bits);
                double p0 = 1.0 - p1;
                double h0 = -(p0*std::log(p0) + p1*std::log(p1))/std::log(2);

                double z = sdsl::size_in_mega_bytes(S) + sdsl::size_in_mega_bytes(D);
                double bz = 8*1024.0*1024.0*z / p_text.size();

                std::cerr << "p0 = " << p0 << std::endl;
                std::cerr << "p1 = " << p1 << std::endl;
                std::cerr << "H0(D) = " << h0 << std::endl;
                std::cerr << "|S| = " << sdsl::size_in_mega_bytes(S) << std::endl;
                std::cerr << "|D| = " << sdsl::size_in_mega_bytes(D) << std::endl;
                std::cerr << "bits/symbol = " << bz << std::endl;
            }

        private:
            size_t m_size;
            std::vector<std::vector<char>> m_codes;
            std::unordered_map<char,std::pair<size_t,size_t>> m_code_index;
            sdsl::bit_vector m_S;
            sdsl::bit_vector m_D;
            sdsl::bit_vector::select_1_type m_select_1_D;
        };
    }
    // namespace seq
}
// namespace varoom

#endif // VAROOM_SEQ_KIMMO_HPP
