#ifndef VAROOM_VCF_VCF_HANDLER_HPP
#define VAROOM_VCF_VCF_HANDLER_HPP

#ifndef VAROOM_UTIL_LAZY_HPP
#include "varoom/util/lazy.hpp"
#endif

#ifndef VAROOM_UTIL_SUBTEXT_HPP
#include "varoom/util/subtext.hpp"
#endif

#ifndef VAROOM_VCF_VCF_INFO_HPP
#include "varoom/vcf/vcf_info.hpp"
#endif

#include <boost/format.hpp>

namespace varoom
{
    namespace vcf
    {
        typedef std::pair<std::string,std::string> key_value_pair;

        class vcf_handler
        {
        public:
            virtual void meta(const std::string& p_key, const std::string& p_value)
            {
            }

            virtual void meta(const std::string& p_kind, const std::vector<key_value_pair>& p_info)
            {
            }

            virtual void operator()(const std::string& p_chr,
                                    const std::int64_t& p_pos,
                                    const std::string& p_id,
                                    const std::string& p_ref,
                                    const std::string& p_alt,
                                    const double& p_qual,
                                    const std::string& p_filter,
                                    const lazy<vcf_info>& p_info,
                                    const lazy<std::vector<vcf_info>>& p_genotypes)
            {
            }

            virtual void error(const size_t& p_line_no, const std::string& p_line, const std::string& p_message)
            {
                throw std::runtime_error(boost::str(boost::format("VCF line %d: %d") % p_line_no % p_message));
            }
        };
    }
    // namespace vcf
}
// namespace varoom

#endif // VAROOM_VCF_VCF_HANDLER_HPP
