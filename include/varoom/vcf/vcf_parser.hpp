#ifndef VAROOM_VCF_VCF_PARSER_HPP
#define VAROOM_VCF_VCF_PARSER_HPP

#ifndef VAROOM_VCF_VCF_HANDLER_HPP
#include "varoom/vcf/vcf_handler.hpp"
#endif

#include <iostream>

namespace varoom
{
    namespace vcf
    {
        class vcf_parser
        {
        public:
            vcf_parser(vcf_handler& p_handler)
                : m_handler(p_handler)
            {
            }

            void parse(std::istream& p_in)
            {
                std::string l;
                std::vector<subtext> parts;
                std::string chr;
                std::int64_t pos;
                std::string id;
                std::string ref;
                std::string alt;
                std::int64_t qual;
                std::string filter;
                std::vector<vcf_info_subtext> genotype_makers;

                size_t line_no = 0;
                while (std::getline(p_in, l))
                {
                    ++line_no;
                    if (starts_with(l, '#'))
                    {
                        continue;
                    }
                    subtext ll(l);
                    parts.clear();
                    ll.split('\t', parts);

                    if (parts.size() < 7)
                    {
                        m_handler.error(line_no, l, "parse error");
                        continue;
                    }
                    chr = parts[0];
                    pos = boost::lexical_cast<std::int64_t>(boost::make_iterator_range(parts[1].first, parts[1].second));
                    id = parts[2];
                    ref = parts[3];
                    alt = parts[4];
                    qual = boost::lexical_cast<double>(boost::make_iterator_range(parts[5].first, parts[5].second));
                    filter = parts[6];
                    vcf_info_subtext info_maker(parts[7]);
                    lazy_vcf_info info([info_maker]() { return info_maker.make(); });

                    genotype_makers.clear();
                    for (size_t i = 9; i < parts.size(); ++i)
                    {
                        genotype_makers.push_back(vcf_info_subtext(parts[8], parts[i]));
                    }
                    lazy<std::vector<vcf_info>> genotypes([genotype_makers]() {
                        std::vector<vcf_info> res;
                        for (size_t i = 0; i < genotype_makers.size(); ++i)
                        {
                            res.push_back(genotype_makers[i].make());
                        }
                        return res;
                    });

                    m_handler(chr, pos, id, ref, alt, qual, filter, info, genotypes);
                }
            }

        private:
            static bool starts_with(const std::string& p_str, char p_ch)
            {
                return p_str.size() > 0 && p_str.front() == p_ch;
            }

            vcf_handler& m_handler;
        };
    }
    // namespace vcf
}
// namespace varoom

#endif // VAROOM_VCF_VCF_PARSER_HPP
