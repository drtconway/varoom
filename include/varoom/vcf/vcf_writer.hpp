#ifndef VAROOM_VCF_VCF_WRITER_HPP
#define VAROOM_VCF_VCF_WRITER_HPP

#ifndef VAROOM_VCF_VCF_HANDLER_HPP
#include "varoom/vcf/vcf_handler.hpp"
#endif

#include <ostream>

namespace varoom
{
    namespace vcf
    {
        class vcf_writer : public vcf_handler
        {
        public:
            vcf_writer(std::ostream& p_out)
                : m_out(p_out)
            {
            }

            virtual void operator()(const std::string& p_chr,
                                    const std::int64_t& p_pos,
                                    const std::string& p_id,
                                    const std::string& p_ref,
                                    const std::string& p_alt,
                                    const std::int64_t& p_qual,
                                    const std::string& p_filter,
                                    const lazy_vcf_info& p_info,
                                    const lazy<std::vector<vcf_info>>& p_genotypes)
            {
                m_out << p_chr
                      << '\t' << p_pos
                      << '\t' << p_id
                      << '\t' << p_ref
                      << '\t' << p_alt
                      << '\t' << p_qual
                      << '\t' << p_filter;
                {
                    const vcf_info& ifo = p_info.get();
                    for (auto itr = ifo.begin(); itr != ifo.end(); ++itr)
                    {
                        m_out << (itr == ifo.begin() ? '\t' : ';');
                        m_out << itr->first << '=' << itr->second;
                    }
                }
                const std::vector<vcf_info>& genotypes = p_genotypes.get();
                if (genotypes.size() > 0)
                {
                    const vcf_info& ifo = genotypes[0];
                    for (auto itr = ifo.begin(); itr != ifo.end(); ++itr)
                    {
                        m_out << (itr == ifo.begin() ? '\t' : ':');
                        m_out << itr->first;
                    }
                }
                for (size_t i = 0; i < genotypes.size(); ++i)
                {
                    const vcf_info& ifo = genotypes[i];
                    for (auto itr = ifo.begin(); itr != ifo.end(); ++itr)
                    {
                        m_out << (itr == ifo.begin() ? '\t' : ':');
                        m_out << itr->second;
                    }
                }
                m_out << std::endl;
            }

            virtual void error(const size_t& p_line_no, const std::string& p_line, const std::string& p_message)
            {
            }

        private:
            std::ostream& m_out;
        };
    }
    // namespace vcf
}
// namespace varoom


#endif // VAROOM_VCF_VCF_WRITER_HPP
