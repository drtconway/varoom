#ifndef VAROOM_HGVS_VCF_TO_HGVSG_HPP
#define VAROOM_HGVS_VCF_TO_HGVSG_HPP

#ifndef VAROOM_HGVS_HGVSG_HANDLER_HPP
#include "varoom/hgvs/hgvsg_handler.hpp"
#endif

#ifndef VAROOM_VCF_VCF_HANDLER_HPP
#include "varoom/vcf/vcf_handler.hpp"
#endif

#ifndef VAROOM_UTIL_SUBTEXT_HPP
#include "varoom/util/subtext.hpp"
#endif

namespace varoom
{
    namespace hgvs
    {
        class vcf_to_hgvsg : public vcf::vcf_handler
        {
        public:
            vcf_to_hgvsg(hgvsg_handler& p_handler)
                : m_handler(p_handler)
            {
            }

            virtual void operator()(const std::string& p_chr,
                                    const std::int64_t& p_pos,
                                    const std::string& p_id,
                                    const std::string& p_ref,
                                    const std::string& p_alt,
                                    const double& p_qual,
                                    const std::string& p_filter,
                                    const varoom::lazy<varoom::vcf::vcf_info>& p_info,
                                    const varoom::lazy<std::vector<varoom::vcf::vcf_info>>& p_genotypes)
            {
                std::vector<subtext> alleles;
                subtext(p_alt).split(',', alleles);
                for (size_t i = 0; i < alleles.size(); ++i)
                {
                    subtext ref = p_ref;
                    subtext alt = alleles[i];
                    if (ref.size() == 1 && alt.size() == 1)
                    {
                        m_handler.sub(p_chr, p_pos, ref, alt);
                        continue;
                    }
                    if (ref.size() == 1 && alt.size() > 1)
                    {
                        ref.first += 1;
                        alt.first += 1;
                        m_handler.ins(p_chr, p_pos, p_pos + 1, alt);
                        continue;
                    }
                    if (ref.size() > 1 && alt.size() == 1)
                    {
                        ref.first += 1;
                        alt.first += 1;
                        m_handler.del(p_chr, p_pos + 1, p_pos + 1 + ref.size());
                        continue;
                    }
                    if (ref.size() > 1 && alt.size() > 1)
                    {
                        int n = 0;
                        while (ref.size() > 0 && alt.size() > 0 && *ref.first == *alt.first)
                        {
                            ref.first += 1;
                            alt.first += 1;
                            n += 1;
                        }
                        if (ref.size() == alt.size() && is_reverse_complement(ref, alt))
                        {
                            m_handler.inv(p_chr, p_pos + n, p_pos + n + ref.size());
                            continue;
                        }
                        m_handler.delins(p_chr, p_pos + n, p_pos + n + ref.size(), alt);
                        continue;
                    }
                    throw std::runtime_error("internal logic error");
                }
            }
        private:
            static bool is_reverse_complement(const subtext& p_lhs, const subtext& p_rhs)
            {
                size_t n = p_lhs.size();
                for (size_t i = 0; i < n; ++i)
                {
                    if (!is_reverse_complement(p_lhs[i], p_rhs[n - i - 1]))
                    {
                        return false;
                    }
                }
                return true;
            }

            static bool is_reverse_complement(char p_lhs, char p_rhs)
            {
                switch (p_lhs)
                {
                    case 'A':
                    case 'a':
                    {
                        return p_rhs == 'T' || p_rhs == 't';
                    }
                    case 'C':
                    case 'c':
                    {
                        return p_rhs == 'G' || p_rhs == 'g';
                    }
                    case 'G':
                    case 'g':
                    {
                        return p_rhs == 'C' || p_rhs == 'c';
                    }
                    case 'T':
                    case 't':
                    {
                        return p_rhs == 'A' || p_rhs == 'a';
                    }
                }
                return false;
            }

            hgvsg_handler& m_handler;
        };
    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_VCF_TO_HGVSG_HPP
