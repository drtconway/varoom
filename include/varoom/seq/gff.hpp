#ifndef VAROOM_SEQ_GFF_HPP
#define VAROOM_SEQ_GFF_HPP

#include <istream>
#include <string>
#include <utility>
#include <vector>
#include <boost/lexical_cast.hpp>

#ifndef VAROOM_UTIL_SUBTEXT_HPP
#include "varoom/util/subtext.hpp"
#endif

namespace varoom
{
    namespace seq
    {
        using attr_and_val = std::pair<std::string,std::string>;

        using handler_function = std::function<void(const std::string& p_seqid, const std::string& p_source, const std::string& p_type,
                                                    const uint64_t& p_start1, const uint64_t& p_end1, const std::string& p_score,
                                                    const std::string& p_strand, const std::string& p_phase,
                                                    const std::vector<attr_and_val>& p_attributes)>;

        class gff3_reader
        {
        public:
            gff3_reader(std::istream& p_in)
                : m_in(p_in)
            {
            }

            template <typename X>
            bool next(X x)
            {
                static_assert(std::is_convertible<X, handler_function>::value,
                              "handler function required.");

                while (true)
                {
                    if (!std::getline(m_in, m_line))
                    {
                        return false;
                    }
                    if (m_line.size() > 0 && m_line[0] == '#')
                    {
                        continue;
                    }
                    subtext st(m_line);
                    st.split('\t', m_parts);
                    if (m_parts.size() != 8 && m_parts.size() != 9)
                    {
                        throw std::runtime_error("wrong number of fields in row");
                    }
                    make<std::string>(m_parts[0], m_seqid);
                    make<std::string>(m_parts[1], m_source);
                    make<std::string>(m_parts[2], m_type);
                    make<uint64_t>(m_parts[3], m_start1);
                    make<uint64_t>(m_parts[4], m_end1);
                    make<std::string>(m_parts[5], m_score);
                    make<std::string>(m_parts[6], m_strand);
                    make<std::string>(m_parts[7], m_phase);
                    m_attributes.clear();
                    if (m_parts.size() == 9)
                    {
                        subtext v = m_parts[8];
                        m_parts.clear();
                        v.split(';', m_parts);
                        attr_and_val w;
                        for (size_t i = 0; i < m_parts.size(); ++i)
                        {
                            subtext u = m_parts[i];
                            m_subparts.clear();
                            u.split('=', m_subparts);
                            make<std::string>(m_subparts[0], w.first);
                            w.second.clear();
                            if (m_subparts.size() > 1)
                            {
                                w.second = std::string(m_subparts[1].first, m_subparts.back().second);
                            }
                            m_attributes.push_back(w);
                        }
                    }
                    x(m_seqid, m_source, m_type, m_start1, m_end1, m_score, m_strand, m_phase, m_attributes);
                    return true;
                }
            }

        private:

            template <typename T>
            void make(const subtext& p_stxt, T& p_res)
            {
                p_res = boost::lexical_cast<T>(boost::make_iterator_range(p_stxt.first, p_stxt.second));
            }

            std::istream& m_in;
            std::string m_line;
            std::vector<subtext> m_parts;
            std::vector<subtext> m_subparts;
            std::string m_seqid;
            std::string m_source;
            std::string m_type;
            uint64_t m_start1;
            uint64_t m_end1;
            std::string m_score;
            std::string m_strand;
            std::string m_phase;
            std::vector<attr_and_val> m_attributes;
        };
    }
    // namespace seq
}
// namespace varoom

#endif // VAROOM_SEQ_GFF_HPP
