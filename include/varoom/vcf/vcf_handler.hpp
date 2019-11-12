#ifndef VAROOM_VCF_VCF_HANDLER_HPP
#define VAROOM_VCF_VCF_HANDLER_HPP

#ifndef VAROOM_UTIL_SUBTEXT_HPP
#include "varoom/util/subtext.hpp"
#endif

#ifndef VAROOM_VCF_VCF_INFO_HPP
#include "varoom/vcf/vcf_info.hpp"
#endif

namespace varoom
{
    namespace vcf
    {
        class vcf_handler
        {
        public:
            virtual void operator()(const std::string& p_chr,
                                    const std::int64_t& p_pos,
                                    const std::string& p_id,
                                    const std::string& p_ref,
                                    const std::string& p_alt,
                                    const std::int64_t& p_qual,
                                    const std::string& p_filter,
                                    const lazy_vcf_info& p_info,
                                    const lazy<std::vector<vcf_info>>& p_genotypes) = 0;

            virtual void error(const size_t& p_line_no, const std::string& p_line, const std::string& p_message) = 0;
        };
    }
    // namespace vcf
}
// namespace varoom

#endif // VAROOM_VCF_VCF_HANDLER_HPP
