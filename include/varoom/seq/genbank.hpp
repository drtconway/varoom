#ifndef VAROOM_SEQ_GENBANK_HPP
#define VAROOM_SEQ_GENBANK_HPP

#include <iostream>
#include <istream>
#include <memory>
#include <regex>
#include <string>
#include <utility>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <nlohmann/json.hpp>

namespace varoom
{
    namespace seq
    {
        typedef nlohmann::json json;

        // Implementation based on
        // https://www.ncbi.nlm.nih.gov/genbank/release/233/

        typedef std::pair<std::string,std::string> string_pair;

        class genbank_location
        {
        public:
            virtual ~genbank_location() {}

            virtual json to_json() const = 0;
        };
        typedef std::shared_ptr<genbank_location> genbank_location_ptr;

        class range_genbank_location : public genbank_location
        {
        public:
            range_genbank_location(const size_t& p_first, const size_t& p_last)
                : m_first(p_first), m_last(p_last)
            {
            }

            const size_t& first() const
            {
                return m_first;
            }

            const size_t& last() const
            {
                return m_last;
            }

            virtual json to_json() const
            {
                json res;
                res.push_back(m_first);
                res.push_back(m_last);
                return res;
            }

        private:
            const size_t m_first;
            const size_t m_last;
        };

        class complement_genbank_location : public genbank_location
        {
        public:
            complement_genbank_location(genbank_location_ptr& p_src)
                : m_src(p_src)
            {
            }

            const genbank_location& source() const
            {
                return *m_src;
            }

            virtual json to_json() const
            {
                json res;
                res["complement"] = source().to_json();
                return res;
            }

        private:
            genbank_location_ptr m_src;
        };

        class join_genbank_location : public genbank_location
        {
        public:
            join_genbank_location(const std::vector<genbank_location_ptr>& p_srcs)
                : m_srcs(p_srcs)
            {
            }

            const std::vector<genbank_location_ptr>& sources() const
            {
                return m_srcs;
            }

            virtual json to_json() const
            {
                json kids;
                for (size_t i = 0; i < m_srcs.size(); ++i)
                {
                    kids.push_back(m_srcs[i]->to_json());
                }
                json res;
                res["join"] = kids;
                return res;
            }

        private:
            std::vector<genbank_location_ptr> m_srcs;
        };

        typedef string_pair genbank_qualifier;

        struct genbank_feature
        {
            std::string type;
            genbank_location_ptr location;
            std::vector<genbank_qualifier> qualifiers;

            json to_json() const
            {
                json res;
                res["type"] = type;
                res["location"] = location->to_json();
                for (size_t i = 0; i < qualifiers.size(); ++i)
                {
                    res["qualifiers"][qualifiers[i].first].push_back(qualifiers[i].second);
                }
                return res;
            }
        };

        struct genbank_entry
        {
            std::string name;
            std::string value;
            std::vector<string_pair> sub_entries;
            std::vector<genbank_feature> features;
            std::string sequence;

            json to_json() const
            {
                json res;
                res[name] = json::array();

                if (value.size() > 0)
                {
                    res[name].push_back(value);
                }

                if (sub_entries.size() > 0)
                {
                    json sub;
                    for (size_t i = 0; i < sub_entries.size(); ++i)
                    {
                        sub[sub_entries[i].first].push_back(sub_entries[i].second);
                    }
                    res[name].push_back(sub);
                }

                if (features.size())
                {
                    json feat;
                    for (size_t i = 0; i < features.size(); ++i)
                    {
                        feat.push_back(features[i].to_json());
                    }
                    res[name].push_back(feat);
                }

                if (sequence.size() > 0)
                {
                    res[name].push_back(sequence);
                }

                return res;
            }
        };

        struct genbank_record
        {
            std::vector<genbank_entry> entries;

            json to_json() const
            {
                json res;
                for (size_t i = 0; i < entries.size(); ++i)
                {
                    res.push_back(entries[i].to_json());
                }
                return res;
            }
        };

        class genbank_reader
        {
        public:
            genbank_reader(std::istream& p_in)
                : m_in(p_in), m_more(true)
            {
                next();
            }

            bool more() const
            {
                return m_more;
            }

            const genbank_record& operator*() const
            {
                return m_record;
            }

            genbank_reader& operator++()
            {
                next();
                return *this;
            }

        private:
            enum line_kind { KW, SUB, FT, NUM, CONT, END };

            void next()
            {
                next_line();
                m_record.entries.clear();
                if (!more())
                {
                    return;
                }
                while (!is_end_line(m_line))
                {
                    parse_entry(m_record.entries);
                }
            }

            void parse_entry(std::vector<genbank_entry>& p_res)
            {
                if (!more() || !is_keyword_line(m_line))
                {
                    if (more())
                    {
                        std::cerr << m_line << std::endl;
                    }
                    throw std::runtime_error("expected KEYWORD");
                }
                string_pair kwl = split_keyword_line(m_line);
                next_line();

                while (more() && is_entry_continuation_line(m_line))
                {
                    boost::algorithm::trim(m_line);
                    kwl.second.push_back('\n');
                    kwl.second.insert(kwl.second.end(), m_line.begin(), m_line.end());
                    next_line();
                }

                std::vector<string_pair> subs;
                while (more() && is_subkeyword_line(m_line))
                {
                    string_pair skwl = split_keyword_line(m_line);
                    next_line();
                    while (more() && is_entry_continuation_line(m_line))
                    {
                        boost::algorithm::trim(m_line);
                        skwl.second.push_back('\n');
                        skwl.second.insert(skwl.second.end(), m_line.begin(), m_line.end());
                        next_line();
                    }
                    subs.push_back(skwl);
                }

                std::vector<genbank_feature> feats;
                while (more() && is_feature_line(m_line))
                {
                    parse_feature(feats);
                }

                std::vector<std::string> seqs;
                while (more() && is_sequence_line(m_line))
                {
                    parse_sequence(seqs);
                }

                genbank_entry e;
                e.name = kwl.first;
                e.value = kwl.second;
                e.sub_entries = subs;
                e.features = feats;
                e.sequence = boost::algorithm::join(seqs, "");
                p_res.push_back(e);
            }

            void parse_feature(std::vector<genbank_feature>& p_features)
            {
                genbank_feature f;
                f.type = std::string(m_line.begin(), m_line.begin() + 20);
                boost::algorithm::trim(f.type);
                m_line.erase(m_line.begin(), m_line.begin() + 20);
                boost::algorithm::trim(m_line);
                f.location = parse_location(m_line);
                next_line();
                while (more() && is_feature_continuation_line(m_line))
                {
                    parse_qualifier(f.qualifiers);
                }
                p_features.push_back(f);
            }

            void parse_qualifier(std::vector<string_pair>& p_qualifiers)
            {
                boost::algorithm::trim(m_line);
                if (!boost::algorithm::starts_with(m_line, "/"))
                {
                    std::cerr << m_line << std::endl;
                    throw std::runtime_error("qualifiers must start with /");
                }
                
                string_pair qual;

                auto p = m_line.begin() + 1;
                auto q = p;
                for (; q != m_line.end(); ++q)
                {
                    if (*q == '=')
                    {
                        break;
                    }
                }
                qual.first = std::string(p, q);
                if (q == m_line.end())
                {
                    p_qualifiers.push_back(qual);
                    next_line();
                    return;
                }
                p = q + 1;
                if (p == m_line.end())
                {
                    // Empty value?
                    next_line();
                    return;
                }
                if (*p != '"')
                {
                    qual.second = std::string(p, m_line.end());
                    p_qualifiers.push_back(qual);
                    next_line();
                    return;
                }
                ++p;
                while (true)
                {
                    while (p == m_line.end())
                    {
                        qual.second.push_back('\n');
                        // string folds on to next line.
                        next_line();
                        if (!more() || !is_feature_continuation_line(m_line))
                        {
                            throw std::runtime_error("unterminated qualifier string value");
                        }
                        m_line.erase(m_line.begin(), m_line.begin() + 21);
                        p = m_line.begin();
                    }
                    if (*p == '"')
                    {
                        ++p;
                        if (p == m_line.end())
                        {
                            p_qualifiers.push_back(qual);
                            next_line();
                            return;
                        }
                        if (*p != '"')
                        {
                            throw std::runtime_error("characters trailing qualifier string value");
                        }
                    }
                    qual.second.push_back(*p);
                    ++p;
                }
            }

            genbank_location_ptr parse_location(const std::string& p_text)
            {
                std::smatch m;

                static const std::regex simple1("([0-9]+)");
                if (std::regex_match(p_text, m, simple1))
                {
                    size_t pos = boost::lexical_cast<size_t>(m[1].str());
                    return genbank_location_ptr(new range_genbank_location(pos, pos));
                }
                static const std::regex simple2("([0-9]+)[.][.]([0-9]+)");
                if (std::regex_match(p_text, m, simple2))
                {
                    size_t first_pos = boost::lexical_cast<size_t>(m[1].str());
                    size_t last_pos = boost::lexical_cast<size_t>(m[2].str());
                    return genbank_location_ptr(new range_genbank_location(first_pos, last_pos));
                }
                static const std::regex uncertain("(<)?([0-9]+)[.][.](>)?([0-9]+)");
                if (std::regex_match(p_text, m, uncertain))
                {
                    std::cerr << "warning: unknown locations not currently supported (" << p_text << ")" << std::endl;
                    size_t first_pos = boost::lexical_cast<size_t>(m[2].str());
                    size_t last_pos = boost::lexical_cast<size_t>(m[4].str());
                    return genbank_location_ptr(new range_genbank_location(first_pos, last_pos));
                }
                std::string u = p_text;
                std::string v;
                int open = 0;
                while (true)
                {
                    for (auto itr = u.begin(); itr != u.end(); ++itr)
                    {
                        switch (*itr)
                        {
                            case '(':
                            {
                                ++open;
                                break;
                            }
                            case ')':
                            {
                                --open;
                                if (open < 0)
                                {
                                    throw std::runtime_error("excess open parentheses");
                                }
                                break;
                            }
                        }
                    }
                    v.insert(v.end(), u.begin(), u.end());
                    if (open == 0)
                    {
                        break;
                    }
                    next_line();
                    if (!more() || !is_feature_continuation_line(m_line))
                    {
                        throw std::runtime_error("unmatched parenthesis in location");
                    }
                    m_line.erase(m_line.begin(), m_line.begin() + 21);
                    u = m_line;
                }
                //std::cerr << v << std::endl;
                static const std::regex complement("complement[(](.*)[)]");
                if (std::regex_match(v, m, complement))
                {
                    std::string inner = m[1].str();
                    genbank_location_ptr locp = parse_location(inner);
                    return genbank_location_ptr(new complement_genbank_location(locp));
                }
                static const std::regex join("join[(](.*)[)]");
                if (std::regex_match(v, m, join))
                {
                    std::vector<genbank_location_ptr> locs;
                    std::string inner = m[1].str();
                    std::regex join_part("([^,]+)");
                    std::sregex_iterator end;
                    for(std::sregex_iterator next(inner.begin(), inner.end(), join_part); next != end; ++next)
                    {
                        m = *next;
                        genbank_location_ptr locp = parse_location(m.str());
                        locs.push_back(locp);
                    }
                    return genbank_location_ptr(new join_genbank_location(locs));
                }
                return genbank_location_ptr();
            }

            void parse_sequence(std::vector<std::string>& p_seqs)
            {
                std::string seq;
                auto p = m_line.begin();
                while (p != m_line.end() && (*p == ' ' || ('0' <= *p && *p <= '9')))
                {
                    ++p;
                }
                for (; p != m_line.end(); ++p)
                {
                    if (*p != ' ' && *p != '\t')
                    {
                        seq.push_back(*p);
                    }
                }
                p_seqs.push_back(seq);
                next_line();
            }

            static bool is_end_line(const std::string& p_line)
            {
                return boost::algorithm::starts_with(p_line, "//");
            }

            static bool is_keyword_line(const std::string& p_line)
            {
                size_t n = count_leading_blanks(p_line);
                return n == 0;
            }

            static bool is_subkeyword_line(const std::string& p_line)
            {
                size_t n = count_leading_blanks(p_line);
                return 1 <= n && n <= 3 && p_line.size() > 12 && p_line[11] == ' ';
            }

            static bool is_feature_line(const std::string& p_line)
            {
                size_t n = count_leading_blanks(p_line);
                return n == 5  && p_line.size() > 20 && p_line[20] == ' ';
            }

            static bool is_sequence_line(const std::string& p_line)
            {
                size_t n = count_leading_blanks(p_line);
                if (n > 9)
                {
                    return false;
                }
                for (auto p = p_line.begin() + n; p != p_line.end() && *p != ' '; ++p)
                {
                    if (!('0' <= *p && *p <= '9'))
                    {
                        return false;
                    }
                }
                return true;
            }

            string_pair split_keyword_line(const std::string& p_line)
            {
                string_pair res;
                auto p = p_line.begin();
                for (; p != p_line.end(); ++p)
                {
                    if (*p != ' ')
                    {
                        break;
                    }
                }
                auto q = p;
                for (; p != p_line.end(); ++p)
                {
                    if (*p == ' ')
                    {
                        break;
                    }
                }
                res.first = std::string(q, p);
                res.second = std::string(p, p_line.end());
                boost::algorithm::trim(res.first);
                boost::algorithm::trim(res.second);
                return res;
            }

            static bool is_entry_continuation_line(const std::string& p_line)
            {
                return count_leading_blanks(p_line) >= 12;
            }

            static bool is_feature_continuation_line(const std::string& p_line)
            {
                return count_leading_blanks(p_line) >= 21;
            }

            static size_t count_leading_blanks(const std::string& p_line)
            {
                size_t n = 0;
                for (auto itr = p_line.begin(); itr != p_line.end() && *itr == ' '; ++itr)
                {
                    n += 1;
                }
                return n;
            }

            static std::string::const_iterator take_while(std::string::const_iterator p_begin, std::string::const_iterator p_end,
                                                          std::function<bool(char)> p_pred)
            {
                for (auto itr = p_begin; itr != p_end; ++itr)
                {
                    if (!p_pred(*itr))
                    {
                        return itr;
                    }
                }
                return p_end;
            }

            void next_line()
            {
                while (more())
                {
                    if (!std::getline(m_in, m_line))
                    {
                        m_more = false;
                    }
                    boost::algorithm::trim_right(m_line);
                    if (m_line.size() > 0)
                    {
                        return;
                    }
                }
            }

            std::istream& m_in;
            bool m_more;
            std::string m_line;
            genbank_record m_record;
        };
    }
    // namespace seq
}
// namespace varoom

#endif // VAROOM_SEQ_GENBANK_HPP
