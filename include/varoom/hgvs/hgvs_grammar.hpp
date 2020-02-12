#ifndef VAROOM_HGVS_HGVS_GRAMMAR_HPP
#define VAROOM_HGVS_HGVS_GRAMMAR_HPP

#include <cctype>
#include <iostream>

#ifndef VAROOM_FUNDS_PARSING_HPP
#include "varoom/funds/parsing.hpp"
#endif

#ifndef VAROOM_HGVS_VARIANT_HPP
#include "varoom/hgvs/variant.hpp"
#endif

namespace varoom
{
    namespace hgvs
    {
        struct accession
        {
            std::string name;
            uint64_t version;
        };

        static varoom::funds::parser<std::string> operator*(varoom::funds::parser<char> p)
        {
            varoom::funds::parser<varoom::funds::list<char>> q = varoom::funds::many0(p);
            return varoom::funds::fmap(varoom::funds::list_to_string, q);
        }

        static varoom::funds::parser<std::string> operator+(varoom::funds::parser<char> p)
        {
            varoom::funds::parser<varoom::funds::list<char>> q = varoom::funds::many1(p);
            return varoom::funds::fmap(varoom::funds::list_to_string, q);
        }

        template <typename T, typename U>
        static varoom::funds::parser<varoom::funds::list<T>> sepby1(varoom::funds::parser<T> p, varoom::funds::parser<U> s)
        {
            using namespace varoom::funds;
            parser<T> q = seq<T>(p, s, [](T t, U u) { return t; });
            return seq<list<T>>(p, many1(q), [](T t, list<T> ts) { return list<T>(t, ts); });
        }

        struct grammar
        {
            /**
             *
             * exprn <- g_var | c_var | n_var | r_var | p_var
             *
             * g_var <- accn ":g." (g_single | g_multi)
             *
             * g_multi <- g_allele (';' g_allele)*
             *
             * g_allele <- '[' g_single (';' g_single)* ']'
             *
             * g_single <- g_id | g_sub | g_ins | g_del | g_delins | g_dup | g_inv | g_con | g_rep
             *
             * g_id <- g_pos '='
             *
             * g_sub <- g_pos nuc '>' nuc
             *
             * g_ins <- g_pos '_' g_pos "ins" (nuc+|num)
             *
             * g_del <- g_loc "del" (nuc+|num)?
             *
             * g_delins <- g_loc "del" (nuc+|num)? "ins" (nuc+|num)
             *
             * g_dup <- g_loc "dup" (nuc+|num)?
             *
             * g_inv <- g_loc "inv" (nuc+|num)?
             *
             * g_con <- g_loc "con" (g_loc | g_ref_loc)
             *
             * g_rep <- g_loc (nuc+ '[' num ']')+
             *
             * g_ref_loc <- (g_ref ':')? g_loc
             *
             * g_ref <- accn
             *
             * g_interval <- g_pos '_' g_pos?
             *
             * g_loc <- g_pos ('_' g_pos)?      == g_interval | g_pos
             *
             * g_pos <- num
             *
             * refseq <- accn ( '(' accn ')' )?
             *
             */

            grammar()
            {
                using namespace varoom::funds;

                using hgvs_parser = parser<variant_ptr>;

            }

            static varoom::funds::parser<variant_ptr> g_var()
            {
                return g_single();
            }

            static varoom::funds::parser<variant_ptr> g_multi()
            {
                using namespace varoom::funds;
                return (seq<genomic_ref>(g_ref(), str(":g."), [](genomic_ref acc, std::string) { return acc; }) >>= [](genomic_ref acc) {
                    return altx({g_phased_alleles(acc), g_unphased_alleles(acc)});
                });
            }

            static varoom::funds::parser<variant_ptr> g_phased_alleles(genomic_ref acc)
            {
                using namespace varoom::funds;
                return seq<variant_ptr>(sepby1(g_phased_allele(acc), sym(';')), [=](list<variant_ptr> aa) {
                    std::vector<variant_ptr> alleles;
                    while (!aa.empty())
                    {
                        alleles.push_back(aa.head());
                        aa = aa.tail();
                    }
                    return variant_ptr(new hgvsg_allele(acc, alleles, true));
                });
            }

            static varoom::funds::parser<variant_ptr> g_unphased_alleles(genomic_ref acc)
            {
                using namespace varoom::funds;
                return seq<variant_ptr>(sepby1(g_change(acc), sym(';')), [=](list<variant_ptr> aa) {
                    std::vector<variant_ptr> alleles;
                    while (!aa.empty())
                    {
                        alleles.push_back(aa.head());
                        aa = aa.tail();
                    }
                    return variant_ptr(new hgvsg_allele(acc, alleles, false));
                });
            }

            static varoom::funds::parser<variant_ptr> g_phased_allele(genomic_ref acc)
            {
                using namespace varoom::funds;
                return seq<variant_ptr>(sym('['), sepby1(g_change(acc), sym(';')), sym(']'), [=](char, list<variant_ptr> aa, char) {
                    std::vector<variant_ptr> alleles;
                    while (!aa.empty())
                    {
                        alleles.push_back(aa.head());
                        aa = aa.tail();
                    }
                    return variant_ptr(new hgvsg_allele(acc, alleles, false));
                });
            }

            static varoom::funds::parser<variant_ptr> g_single()
            {
                using namespace varoom::funds;
                return (seq<genomic_ref>(g_ref(), str(":g."),
                    [](genomic_ref acc, std::string) { return acc; }) >>= [](genomic_ref acc) {
                        return g_change(acc);
                    });
            }

            static varoom::funds::parser<variant_ptr> g_change(genomic_ref acc)
            {
                using namespace varoom::funds;
                return altx({
                    g_id(acc), g_sub(acc),
                    g_delins(acc),
                    g_del(acc), g_ins(acc)
                });
            }

            static varoom::funds::parser<variant_ptr> g_id(genomic_ref acc)
            {
                using namespace varoom::funds;
                return seq<variant_ptr>(g_pos(), sym('='), [=](genomic_pos pos, char) {
                    return variant_ptr(new hgvsg_id(acc, pos));
                });
            }

            static varoom::funds::parser<variant_ptr> g_sub(genomic_ref acc)
            {
                using namespace varoom::funds;
                return seq<variant_ptr>(g_pos(), nuc(), sym('>'), nuc(), [=](genomic_pos pos, nucleotide ref, char, nucleotide alt) {
                    return variant_ptr(new hgvsg_sub(acc, pos, ref, alt));
                });
            }

            static varoom::funds::parser<variant_ptr> g_ins(genomic_ref acc)
            {
                using namespace varoom::funds;
                return seq<variant_ptr>(g_interval(), str("ins"), +nuc(), [=](genomic_locus loc, std::string, nucleotides alt) {
                    return variant_ptr(new hgvsg_ins(acc, loc, alt));
                });
            }

            static varoom::funds::parser<variant_ptr> g_del(genomic_ref acc)
            {
                using namespace varoom::funds;
                return seq<variant_ptr>(g_loc(), str("del"), *nuc(), [=](genomic_locus loc, std::string, nucleotides ref) {
                    return variant_ptr(new hgvsg_del(acc, loc, ref));
                });
            }

            static varoom::funds::parser<variant_ptr> g_delins(genomic_ref acc)
            {
                using namespace varoom::funds;
                return seq<variant_ptr>(g_loc(), str("del"), *nuc(), str("ins"), +nuc(), [=](genomic_locus loc, std::string, nucleotides ref, std::string, nucleotides alt) {
                    return variant_ptr(new hgvsg_delins(acc, loc, ref, alt));
                });
            }

            static varoom::funds::parser<variant_ptr> g_dup(genomic_ref acc)
            {
                using namespace varoom::funds;
                return seq<variant_ptr>(g_loc(), str("dup"), *nuc(), [=](genomic_locus loc, std::string s, nucleotides ref) {
                    return variant_ptr(new hgvsg_dup(acc, loc, ref));
                });
            }

            static varoom::funds::parser<variant_ptr> g_inv(genomic_ref acc)
            {
                using namespace varoom::funds;
                return seq<variant_ptr>(g_loc(), str("inv"), *nuc(), [=](genomic_locus loc, std::string s, nucleotides ref) {
                    return variant_ptr(new hgvsg_inv(acc, loc, ref));
                });
            }

            using ref_and_loc = std::pair<genomic_ref,genomic_locus>;
            static varoom::funds::parser<variant_ptr> g_con(genomic_ref acc)
            {
                using namespace varoom::funds;
                return seq<variant_ptr>(g_loc(), str("con"), g_ref_loc(acc),
                    [=](genomic_locus loc, std::string s, ref_and_loc oth) {
                        return variant_ptr(new hgvsg_con(acc, loc, oth.first, oth.second));
                });
            }

            static varoom::funds::parser<genomic_locus> g_loc()
            {
                using namespace varoom::funds;
                return altx({
                    g_interval(),
                    fmap([](genomic_pos p) {
                        return genomic_locus(p, p);
                    }, g_pos())
                });
            }

            static varoom::funds::parser<genomic_locus> g_interval()
            {
                using namespace varoom::funds;
                return seq<genomic_locus>(g_pos(), sym('_'), g_pos(), [](genomic_pos f, char c, genomic_pos l) {
                    return genomic_locus(f, l);
                });
            }

            static varoom::funds::parser<genomic_pos> g_pos()
            {
                using namespace varoom::funds;
                return fmap([](int64_t x) { return genomic_pos(uint64_t(x)); }, num());
            }

            static varoom::funds::parser<ref_and_loc> g_ref_loc(genomic_ref acc)
            {
                using namespace varoom::funds;
                return altx({
                    seq<ref_and_loc>(g_ref(), sym(':'), g_loc(), [](genomic_ref ref, char, genomic_locus loc) {
                        return ref_and_loc(ref, loc);
                    }),
                    seq<ref_and_loc>(g_loc(), [=](genomic_locus loc) {
                        return ref_and_loc(acc, loc);
                    })
                });
            }

            static varoom::funds::parser<genomic_ref> g_ref()
            {
                using namespace varoom::funds;
                return fmap([](accession x) {
                    return genomic_ref{x.name, x.version};
                }, accn());
            }

            static varoom::funds::parser<accession> accn()
            {
                using namespace varoom::funds;
                std::function<std::string(char,std::string)> f = &cat;
                parser<std::string> part1 = seq<std::string>(alpha(), *alnum(), f);
                parser<std::string> part2 = seq<std::string>(sym('_'), *alnum(), f);
                parser<uint64_t> part3 = seq<uint64_t>(sym('.'), num(), [](char c, uint64_t n) {
                    return n;
                });
                return altx({
                    seq<accession>(part1, part2, part3, [](std::string p1, std::string p2, uint64_t v) {
                        return accession{p1+p2, v};
                    }),
                    seq<accession>(part1, part2, [](std::string p1, std::string p2) {
                        return accession{p1+p2, 0};
                    }),
                    seq<accession>(part1, [](std::string p1) {
                        return accession{p1, 0};
                    }),
                });
            }

            static varoom::funds::parser<nucleotide> nuc()
            {
                using namespace varoom::funds;
                parser<char> s = sym();
                return sat(s, [](char c) {
                    switch (c)
                    {
                        case 'A':
                        case 'a':
                        case 'C':
                        case 'c':
                        case 'G':
                        case 'g':
                        case 'T':
                        case 't':
                        case 'N':
                        case 'n':
                            return true;
                        default:
                            return false;
                    }
                });
            }

            static varoom::funds::parser<char> alpha()
            {
                using namespace varoom::funds;
                parser<char> s = sym();
                return sat(s, [](char c) { return std::isalpha(c); });
            }

            static varoom::funds::parser<char> alnum()
            {
                using namespace varoom::funds;
                parser<char> s = sym();
                return sat(s, [](char c) { return std::isalnum(c); });
            }

            static varoom::funds::parser<char> digit()
            {
                using namespace varoom::funds;
                parser<char> s = sym();
                return sat(s, [](char c) { return std::isdigit(c); });
            }

            static varoom::funds::parser<uint64_t> num()
            {
                using namespace varoom::funds;
                return fmap([](std::string ds) {
                    uint64_t x = 0;
                    for (auto d = ds.begin(); d != ds.end(); ++d)
                    {
                        x = (10 * x) + uint64_t(*d - '0');
                    }
                    return x;
                }, +digit());
            }

            static std::string cat(char c, std::string s)
            {
                std::string r;
                r.reserve(1+s.size());
                r.push_back(c);
                r.insert(r.end(), s.begin(), s.end());
                return r;
            }
        };
    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_HGVS_GRAMMAR_HPP
