#ifndef VAROOM_VCF_VCF_HANDLER_HPP
#define VAROOM_VCF_VCF_HANDLER_HPP

#ifndef VAROOM_UTIL_SUBTEXT_HPP
#include "varoom/util/subtext.hpp"
#endif

namespace varoom
{
    namespace vcf
    {
        class vcf_info
        {
        public:
            vcf_info(const subtext& p_source)
            {
                m_sources.push_back(p_source);
            }

            vcf_info(const subtext& p_key_source, p_val_source)
            {
                m_sources.push_back(p_key_source);
                m_sources.push_back(p_val_source);
            }

            size_t size() const
            {
                scan_if_necessary();
            }

            void keys(std::vector<std::string>& p_keys) const
            {
                scan_if_necessary();
                p_keys.insert(p_keys.end(), m_keys.begin(), m_keys.end());
            }

            template <typename T>
            bool get(const std::string& p_id, T& p_val) const
            {
                scan_if_necessary();
                size_t i = find(p_id);
                if (i < 0)
                {
                    return false;
                }
                p_val = boost::lexical_cast<T>(static_cast<std::string>(m_vals[i]));
                return true;
            }

            template <typename T>
            bool get(const std::string& p_id, std::vector<T>& p_vals, char p_sep) const
            {
                scan_if_necessary();
                size_t i = find(p_id);
                if (i < 0)
                {
                    return false;
                }
                std::vector<subtext> parts;
                m_vals.split(p_sep, parts);
                for (size_t j = 0; j < parts.size(); ++j)
                {
                    p_vals.push_back(boost::lexical_cast<T>(static_cast<std::string>(parts[j])));
                }
                return true;
            }

        private:
            size_t find(const std::string& p_id) const
            {
                scan_if_necessary();
                for (size_t i = 0; i < m_keys.size(); ++i)
                {
                    if (m_keys[i] == p_id)
                    {
                        return i;
                    }
                }
                return -1;
            }

            void scan_if_necessary() const
            {
                if (m_keys.size() > 0)
                {
                    return;
                }
                if (m_sources.size() == 1)
                {
                    std::vector<subtext> parts;
                    std::vector<subtext> subparts;

                    m_source.split(';', parts);
                    for (size_t i = 0; i < parts.size(); ++i)
                    {
                        subparts.clear();
                        parts[i].split('=', subparts);
                        if (subparts.size() < 2)
                        {
                            // malformed.
                            // TODO: figure out a better strategy than just skipping it.
                            //
                            continue;
                        }
                        if (subparts.size() > 2)
                        {
                            subparts[1].second = subparts.back().second;
                        }
                        m_keys.push_back(subparts[0]);
                        m_vals.push_back(subparts[1]);
                    }
                }
                else
                {
                }
            }

            subtext m_source;
            mutable std::vector<subtext> m_keys;
            mutable std::vector<subtext> m_vals;
        };

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
                                    const vcf_info& p_info,
                                    const std::vector<vcf_info>& p_genotypes) = 0;

            virtual void error(const size_t& p_line_no, const std::string& p_line, const std::string& p_message) = 0;
        };
    }
    // namespace vcf
}
// namespace varoom

#endif // VAROOM_VCF_VCF_HANDLER_HPP
