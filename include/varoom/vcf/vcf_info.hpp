#ifndef VAROOM_VCF_VCF_INFO_HPP
#define VAROOM_VCF_VCF_INFO_HPP

#include <boost/lexical_cast.hpp>
#include <boost/range/iterator_range.hpp>

#ifndef VAROOM_UTIL_SUBTEXT_HPP
#include "varoom/util/subtext.hpp"
#endif

#ifndef VAROOM_UTIL_LAZY_HPP
#include "varoom/util/lazy.hpp"
#endif

#ifndef VAROOM_UTIL_ORDERED_ASSOCIATION_HPP
#include "varoom/util/ordered_association.hpp"
#endif

namespace varoom
{
    namespace vcf
    {
        using vcf_info = ordered_association;

        class vcf_info_subtext
        {
        public:
            vcf_info_subtext()
            {
            }

            vcf_info_subtext(const subtext& p_source)
            {
                m_sources.push_back(p_source);
            }

            vcf_info_subtext(const subtext& p_key_source, const subtext& p_val_source)
            {
                m_sources.push_back(p_key_source);
                m_sources.push_back(p_val_source);
            }

            vcf_info make() const
            {
                ordered_association res;

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
                        res[subparts[0]] = subparts[1];
                    }
                }
                else
                {
                    std::vector<subtext> key_parts;
                    std::vector<subtext> val_parts;

                    m_sources[0].split(':', key_parts);
                    m_sources[1].split(':', val_parts);
                    
                    if (key_parts.size() != val_parts.size())
                    {
                        throw std::runtime_error("format & genotype fields do not match");
                    }

                    for (size_t i = 0; i < key_parts.size(); ++i)
                    {
                        res[key_parts[i]] = val_parts[i];
                    }
                }

                return res;
            }

        private:
            std::vector<subtext> m_sources;
        };
    }
    // namespace vcf
}
// namespace varoom

#endif // VAROOM_VCF_VCF_INFO_HPP
