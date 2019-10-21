#ifndef VAROOM_HGVS_HGVSG_PARSER_HPP
#define VAROOM_HGVS_HGVSG_PARSER_HPP

#ifndef VAROOM_HGVS_HGVSG_HANDLER_HPP
#include "varoom/hgvs/hgvsg_handler.hpp"
#endif

#include <initializer_list>
#include <regex>
#include <boost/lexical_cast.hpp>

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
                //const std::regex pieces_regex("([a-z]+)\\.([a-z]+)");
                //std::smatch pieces_match;
                //if (std::regex_match(fname, pieces_match, pieces_regex)) {
                //    std::cout << fname << '\n';
                //    for (size_t i = 0; i < pieces_match.size(); ++i) {
                //        std::ssub_match sub_match = pieces_match[i];
                //        std::string piece = sub_match.str();
                //        std::cout << "  submatch " << i << ": " << piece << '\n';
                //    }   
                //}   
            }
        private:
            static constexpr const char* accPat = "([^:]+)";
            static constexpr const char* gPosPat = "([0-9]+)";
            static constexpr const char* cPosPat = "(([-*]?)([0-9]+)([-+]?[0-9]+)?)";
            static constexpr const char* pPosPat = "([0-9]+)";
            static constexpr const char* aaPat = "([A-Z][a-z][a-z]|[*])";

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
            //static const char* hgvsgDel1 = Pattern.compile(accPat + ":g[.]" + gPosPat + "del" + "([ACGTacgtNn]*)");
            //static const char* hgvsgDel2 = Pattern.compile(accPat + ":g[.]" + gPosPat + "_" + gPosPat + "del" + "([ACGTacgtNn]*)");
            //static const char* hgvsgDelIns1 = Pattern.compile(accPat + ":g[.]" + gPosPat + "delins" + "([ACGTacgt]+)");
            //static const char* hgvsgDelIns2 = Pattern.compile(accPat + ":g[.]" + gPosPat + "_" + gPosPat + "delins" + "([ACGTacgt]+)");
            //static const char* hgvsgDup1 = Pattern.compile(accPat + ":g[.]" + gPosPat + "dup" + "([ACGTacgt]*)");
            //static const char* hgvsgDup2 = Pattern.compile(accPat + ":g[.]" + gPosPat + "_" + gPosPat + "dup" + "([ACGTacgt]*)");
            //static const char* hgvsgInv = Pattern.compile(accPat + ":g[.]" + gPosPat + "_" + gPosPat + "inv" + "([ACGTacgt]*)");
            //static const char* hgvsgRep = Pattern.compile(accPat + ":g[.]" + gPosPat + "([ACGTacgt]+)\\[([0-9]+)\\]");
            //static const char* hgvsgSil1 = Pattern.compile(accPat + ":g[.]" + gPosPat + "=" + "([ACGTacgt]?)");
            //static const char* hgvsgSil2 = Pattern.compile(accPat + ":g[.]" + gPosPat + "_" + gPosPat + "=" + "([ACGTacgt]*)");
        };
    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_HGVSG_PARSER_HPP
