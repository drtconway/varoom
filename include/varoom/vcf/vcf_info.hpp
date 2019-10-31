#ifndef VAROOM_VCF_VCF_INFO_HPP
#define VAROOM_VCF_VCF_INFO_HPP

#include <boost/lexical_cast.hpp>

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
            vcf_info()
            {
            }

            vcf_info(const subtext& p_source)
            {
                m_sources.push_back(p_source);
            }

            vcf_info(const subtext& p_key_source, const subtext& p_val_source)
            {
                m_sources.push_back(p_key_source);
                m_sources.push_back(p_val_source);
            }

            size_t size() const
            {
                scan_if_necessary();
                return m_keys.size();
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
                int i = find(p_id);
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
                int i = find(p_id);
                if (i < 0)
                {
                    return false;
                }
                std::vector<subtext> parts;
                m_vals[i].split(p_sep, parts);
                for (size_t j = 0; j < parts.size(); ++j)
                {
                    p_vals.push_back(boost::lexical_cast<T>(static_cast<std::string>(parts[j])));
                }
                return true;
            }

        private:
            int find(const std::string& p_id) const
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

                    m_sources[0].split(';', parts);
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
                    std::vector<subtext> key_parts;
                    std::vector<subtext> val_parts;
                    std::vector<subtext> subparts;

                    m_sources[0].split(':', key_parts);
                    m_sources[1].split(':', val_parts);
                    
                    for (size_t i = 0; i < key_parts.size(); ++i)
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
            }

            std::vector<subtext> m_sources;
            mutable std::vector<subtext> m_keys;
            mutable std::vector<subtext> m_vals;
        };
    }
    // namespace vcf
}
// namespace varoom

#endif // VAROOM_VCF_VCF_INFO_HPP
