#ifndef VAROOM_SEQ_GENBANK_HPP
#define VAROOM_SEQ_GENBANK_HPP

#include <istream>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <boost/algorithm/string.hpp>

namespace varoom
{
    namespace seq
    {
        typedef std::pair<std::string,std::string> string_pair;

        class genbank_location
        {
        public:
            virtual ~genbank_location() {}
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

        private:
            std::vector<genbank_location_ptr> m_srcs;
        };

        typedef string_pair genbank_qualifier;

        struct genbank_feature
        {
            std::string type;
            genbank_location_ptr location;
            std::vector<genbank_qualifier> qualifiers;
        };

        struct genbank_entry
        {
            std::string name;
            std::string value;
            std::vector<string_pair> sub_entries;
            std::vector<genbank_feature> features;
        };

        struct genbank_record
        {
            std::vector<genbank_entry> entries;
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
                    throw std::runtime_error("expected KEYWORD");
                }
                string_pair kwl = split_keyword_line(m_line);
                next_line();
                while (more() && is_continuation_line(m_line))
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
                    while (more() && is_continuation_line(m_line))
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
                genbank_entry e;
                e.name = kwl.first;
                e.value = kwl.second;
                e.sub_entries = subs;
                e.features = feats;
                p_res.push_back(e);
            }

            void parse_feature(std::vector<genbank_feature>& p_features)
            {
                throw std::runtime_error("not implemented");
            }

            static bool is_end_line(const std::string& p_line)
            {
                return boost::algorithm::starts_with(p_line, "//");
            }

            static bool is_keyword_line(const std::string& p_line)
            {
                size_t n = count_leading_blanks(p_line);
                return n == 0 && p_line.size() > 12 && p_line[11] == ' ';
            }

            static bool is_subkeyword_line(const std::string& p_line)
            {
                size_t n = count_leading_blanks(p_line);
                return 1 <= n && n <= 3 && p_line.size() > 12 && p_line[11] == ' ';
            }

            static bool is_feature_line(const std::string& p_line)
            {
                size_t n = count_leading_blanks(p_line);
                return n == 5  && p_line.size() > 12 && p_line[11] == ' ';
            }

            string_pair split_keyword_line(const std::string& p_line)
            {
                string_pair res;
                res.first = std::string(p_line.begin(), p_line.begin() + 12);
                res.second = std::string(p_line.begin() + 12, p_line.end());
                boost::algorithm::trim(res.first);
                boost::algorithm::trim(res.second);
                return res;
            }

            static bool is_continuation_line(const std::string& p_line)
            {
                return count_leading_blanks(p_line) == 12;
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

            void next_line()
            {
                if (!std::getline(m_in, m_line))
                {
                    m_more = false;
                }
                boost::algorithm::trim_right(m_line);
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
