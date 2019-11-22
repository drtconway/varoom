#ifndef VAROOM_VCF_VCF_PARSER_HPP
#define VAROOM_VCF_VCF_PARSER_HPP

#ifndef VAROOM_VCF_VCF_HANDLER_HPP
#include "varoom/vcf/vcf_handler.hpp"
#endif

#ifndef VAROOM_UTIL_TEXT_HPP
#include "varoom/util/text.hpp"
#endif

#include <iostream>

namespace varoom
{
    namespace vcf
    {
        namespace detail
        {
            enum tok_kind { WORD, EQ, COMMA, STR };
            struct tok
            {
                tok_kind kind;
                std::string value;

                tok(tok_kind p_kind)
                    : kind(p_kind)
                {
                }

                tok(tok_kind p_kind, const std::string& p_value)
                    : kind(p_kind), value(p_value)
                {
                }

                static bool scanner(const subtext& p_txt, std::vector<tok>& p_toks)
                {
                    auto cur = p_txt.first;
                    while (cur != p_txt.second)
                    {
                        char c = *cur;
                        switch (c)
                        {
                            case '=':
                            {
                                p_toks.push_back(tok(detail::EQ));
                                ++cur;
                                break;
                            }
                            case ',':
                            {
                                p_toks.push_back(tok(detail::COMMA));
                                ++cur;
                                break;
                            }
                            case '"':
                            {
                                std::string str;
                                ++cur;
                                while (cur != p_txt.second && *cur != '"')
                                {
                                    if (*cur == '\\')
                                    {
                                        ++cur;
                                        if (cur == p_txt.second)
                                        {
                                            return false;
                                        }
                                    }
                                    str.push_back(*cur);
                                    ++cur;
                                }
                                if (cur == p_txt.second)
                                {
                                    return false;
                                }
                                ++cur;
                                p_toks.push_back(tok(detail::STR, str));
                                break;
                            }
                            default:
                            {
                                std::string val;
                                while (cur != p_txt.second && *cur != '"' && *cur != ',' && *cur != '=')
                                {
                                    val.push_back(*cur);
                                    ++cur;
                                }
                                p_toks.push_back(tok(detail::WORD, val));
                                break;
                            }
                        }
                    }
                    return true;
                }
            };
        }
        // namespace detail

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
                double qual;
                std::string filter;
                std::vector<vcf_info_subtext> genotype_makers;

                size_t line_no = 0;
                while (std::getline(p_in, l))
                {
                    ++line_no;
                    if (text::starts_with(l, "##"))
                    {
                        parse_meta(l);
                        continue;
                    }
                    if (text::starts_with(l, '#'))
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
                    lazy<vcf_info> info([info_maker]() { return info_maker.make(); });

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
            void parse_meta(const std::string& p_line)
            {
                std::vector<subtext> parts;
                // Skip over the ##
                subtext ll(p_line.begin() + 2, p_line.end());
                ll.split('=', parts);

                if (parts.size() < 2)
                {
                    return;
                }

                subtext kind = parts[0];
                subtext value = subtext(parts[1].first, parts.back().second);

                if (text::starts_with(value, '<') && text::ends_with(value, '>'))
                {
                    value.first += 1;
                    value.second -= 1;
                    std::vector<key_value_pair> kvs;
                    if (!parse_key_value_pairs(value, kvs))
                    {
                        m_handler.meta(kind, value);
                    }
                    m_handler.meta(kind, kvs);
                }
                else
                {
                    m_handler.meta(kind, value);
                }
            }

            static bool parse_key_value_pairs(const subtext& p_str, std::vector<key_value_pair>& p_kvs)
            {
                std::vector<detail::tok> toks;
                if (!detail::tok::scanner(p_str, toks))
                {
                    return false;
                }
                size_t i = 0;
                while (i < toks.size())
                {
                    key_value_pair kv;
                    if (toks[i].kind != detail::WORD)
                    {
                        return false;
                    }
                    kv.first = toks[i].value;
                    ++i;
                    if (i == toks.size())
                    {
                        return false;
                    }
                    if (toks[i].kind == detail::COMMA)
                    {
                        p_kvs.push_back(kv);
                        ++i;
                        continue;
                    }
                    if (toks[i].kind != detail::EQ)
                    {
                        return false;
                    }
                    ++i;
                    if (i == toks.size() || (toks[i].kind != detail::WORD && toks[i].kind != detail::STR))
                    {
                        return false;
                    }
                    kv.second = toks[i].value;
                    p_kvs.push_back(kv);
                    ++i;
                    if (i == toks.size())
                    {
                        continue;
                    }
                    if (toks[i].kind != detail::COMMA)
                    {
                        return false;
                    }
                    ++i;
                }
                return true;
            }

            vcf_handler& m_handler;
        };
    }
    // namespace vcf
}
// namespace varoom

#endif // VAROOM_VCF_VCF_PARSER_HPP
