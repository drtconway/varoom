#ifndef VAROOM_VCF_VCF_PARSER_HPP
#define VAROOM_VCF_VCF_PARSER_HPP

#ifndef VAROOM_VCF_VCF_HANDLER_HPP
#include "varoom/vcf/vcf_handler.hpp"
#endif

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
                vector<subtext> parts;
                std::string chr;
                std::int64_t pos;
                std::string id;
                std::string ref;
                std::string alt;
                std::int64_t qual;
                std::string filter;
                vcf_info info;
                std::vector<vcf_info> genotypes;

                size_t line_no = 0;
                while (std::getline(p_in, l))
                {
                    ++line_no;
                    if (l.starts_with('#'))
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
                    if (!boost::conversion::type_lexical_convert(parts[1], pos))
                    {
                        m_handler.error(line_no, l, "could not convert position");
                        continue;
                    }
                    id = parts[2];
                    ref = parts[3];
                    alt = parts[4];
                    qual = parts[5];
                    filter = parts[6];
                    info = vcf_info(parts[7]);

                    genotypes.clear();
                    for (size_t i = 9; i < parts.size(); ++i)
                    {
                        genotypes.push_back(vcf_info(parts[8], parts[i]));
                    }

                    m_handler(chr, pos, id, ref, alt, qual, filter, info, genotypes);
                }
            }

        private:
            vcf_handler& m_handler;
        };
    }
    // namespace vcf
}
// namespace varoom

#endif // VAROOM_VCF_VCF_PARSER_HPP
