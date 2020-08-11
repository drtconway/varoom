#ifndef VAROOM_HGVS_HGVSG_PARSER_HPP
#define VAROOM_HGVS_HGVSG_PARSER_HPP

#ifndef VAROOM_HGVS_HGVSG_HANDLER_HPP
#include "varoom/hgvs/hgvsg_handler.hpp"
#endif

#include <initializer_list>
#include <regex>
#include <boost/lexical_cast.hpp>
#include <boost/range/iterator_range.hpp>

namespace varoom
{
    namespace hgvs
    {
        class hgvsg_parser
        {
        public:
            static void parse(const std::string& p_txt, hgvsg_handler& p_handler)
            {
                std::smatch m;
                if (std::regex_match(p_txt, m, hgvsgSub()))
                {
                    p_handler.sub(m[1], boost::lexical_cast<std::int64_t>(m[2]), m[3], m[4]);
                    return;
                }
                if (std::regex_match(p_txt, m, hgvsgIns()))
                {
                    p_handler.ins(m[1], boost::lexical_cast<std::int64_t>(m[2]), boost::lexical_cast<std::int64_t>(m[3]), m[4]);
                    return;
                }
                if (std::regex_match(p_txt, m, hgvsgDel1()))
                {
                    std::int64_t first = boost::lexical_cast<std::int64_t>(m[2]);
                    p_handler.del(m[1], first, first);
                    return;
                }
                if (std::regex_match(p_txt, m, hgvsgDel2()))
                {
                    std::int64_t first = boost::lexical_cast<std::int64_t>(m[2]);
                    std::int64_t last = boost::lexical_cast<std::int64_t>(m[3]);
                    p_handler.del(m[1], first, last);
                    return;
                }
                if (std::regex_match(p_txt, m, hgvsgDelIns1()))
                {
                    std::int64_t first = boost::lexical_cast<std::int64_t>(m[2]);
                    p_handler.delins(m[1], first, first, m[3]);
                    return;
                }
                if (std::regex_match(p_txt, m, hgvsgDelIns2()))
                {
                    std::int64_t first = boost::lexical_cast<std::int64_t>(m[2]);
                    std::int64_t last = boost::lexical_cast<std::int64_t>(m[3]);
                    p_handler.delins(m[1], first, last, m[4]);
                    return;
                }
                if (std::regex_match(p_txt, m, hgvsgDup1()))
                {
                    std::int64_t first = boost::lexical_cast<std::int64_t>(m[2]);
                    p_handler.dup(m[1], first, first);
                    return;
                }
                if (std::regex_match(p_txt, m, hgvsgDup2()))
                {
                    std::int64_t first = boost::lexical_cast<std::int64_t>(m[2]);
                    std::int64_t last = boost::lexical_cast<std::int64_t>(m[3]);
                    p_handler.dup(m[1], first, last);
                    return;
                }
                if (std::regex_match(p_txt, m, hgvsgInv()))
                {
                    std::int64_t first = boost::lexical_cast<std::int64_t>(m[2]);
                    std::int64_t last = boost::lexical_cast<std::int64_t>(m[3]);
                    p_handler.inv(m[1], first, last);
                    return;
                }
            }
        private:

            static constexpr const char* accPat = "([^:]+)";
            static constexpr const char* gPosPat = "([0-9]+)";

            static std::string compose(std::initializer_list<const char*> p_parts)
            {
                std::string s;
                for (auto i = p_parts.begin(); i != p_parts.end(); ++i)
                {
                    s.insert(s.size(), *i);
                }
                return s;
            }

            static const std::regex& hgvsgSub()
            {
                //static std::regex r(std::string(accPat) + std::string(":g[.]") + std::string(gPosPat) + std::string("([ACGTacgt])>([ACGTacgt])"));
                static std::regex r(compose({accPat, ":g[.]", gPosPat, "([ACGTacgt])>([ACGTacgt])"}));
                return r;
            }

            static const std::regex& hgvsgIns()
            {
                static std::regex r(compose({accPat, ":g[.]", gPosPat, "_", gPosPat, "ins", "([ACGTacgt]+)"}));
                return r;
            }

            static const std::regex& hgvsgDel1()
            {
                static std::regex r(compose({accPat, ":g[.]", gPosPat, "del", "([ACGTacgtNn]*)"}));
                return r;
            }
            static const std::regex& hgvsgDel2()
            {
                static std::regex r(compose({accPat, ":g[.]", gPosPat, "_", gPosPat, "del", "([ACGTacgtNn]*)"}));
                return r;
            }

            static const std::regex& hgvsgDelIns1()
            {
                static std::regex r(compose({accPat, ":g[.]", gPosPat, "delins", "([ACGTacgt]+)"}));
                return r;
            }
            static const std::regex& hgvsgDelIns2()
            {
                static std::regex r(compose({accPat, ":g[.]", gPosPat, "_", gPosPat, "delins", "([ACGTacgt]+)"}));
                return r;
            }

            static const std::regex& hgvsgDup1()
            {
                static std::regex r(compose({accPat, ":g[.]", gPosPat, "dup", "([ACGTacgt]*)"}));
                return r;
            }
            static const std::regex& hgvsgDup2()
            {
                static std::regex r(compose({accPat, ":g[.]", gPosPat, "_", gPosPat, "dup", "([ACGTacgt]*)"}));
                return r;
            }

            static const std::regex& hgvsgInv()
            {
                static std::regex r(compose({accPat, ":g[.]", gPosPat, "_", gPosPat, "inv", "([ACGTacgt]*)"}));
                return r;
            }

            //static const char* hgvsgRep = Pattern.compile(accPat + ":g[.]" + gPosPat + "([ACGTacgt]+)\\[([0-9]+)\\]");
            //static const char* hgvsgSil1 = Pattern.compile(accPat + ":g[.]" + gPosPat + "=" + "([ACGTacgt]?)");
            //static const char* hgvsgSil2 = Pattern.compile(accPat + ":g[.]" + gPosPat + "_" + gPosPat + "=" + "([ACGTacgt]*)");
        };
    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_HGVSG_PARSER_HPP
