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

        struct grammar_basics
        {
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

            static varoom::funds::parser<spliced_position> spos()
            {
                using namespace varoom::funds;
                parser<int64_t> utr5 = seq<int64_t>(sym('-'), num(), [](char, uint64_t x) { return -int64_t(x); });
                parser<int64_t> utr3 = seq<int64_t>(sym('*'), num(), [](char, uint64_t x) { return int64_t(x); });
                parser<int64_t> off = seq<int64_t>(oneof("+-"), num(), [](char s, uint64_t x) {
                    int64_t y = x;
                    return (s == '+' ? y : -y);
                });

                return seq<spliced_position>(altx({
                        fmap([](int64_t x) { return spliced_position{UTR5, x, 0}; }, utr5),
                        fmap([](int64_t x) { return spliced_position{UTR3, x, 0}; }, utr3),
                        fmap([](uint64_t x) { return spliced_position{TR, int64_t(x), 0}; }, num())
                    }), optional(off), [](spliced_position r, maybe<int64_t> moff) {
                        if (!moff.nothing())
                        {
                            r.off = moff.just();
                        }
                        return r;
                    });
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

            static std::string cat(char c, std::string s)
            {
                std::string r;
                r.reserve(1+s.size());
                r.push_back(c);
                r.insert(r.end(), s.begin(), s.end());
                return r;
            }
        };

        template <typename X>
        struct grammar_traits {};

        template <>
        struct grammar_traits<hgvsg>
        {
            static varoom::funds::parser<hgvsg::position> pos()
            {
                using namespace varoom::funds;
                return fmap([](int64_t x) { return hgvsg::position(uint64_t(x)); }, grammar_basics::num());
            }

            static varoom::funds::parser<hgvsg::reference> ref()
            {
                using namespace varoom::funds;
                return fmap([](accession a) { return hgvsg::reference{a.name, a.version}; }, grammar_basics::accn());
            }
        };

        template <>
        struct grammar_traits<hgvsc>
        {
            static varoom::funds::parser<hgvsc::position> pos()
            {
                using namespace varoom::funds;
                return grammar_basics::spos();
            }

            static varoom::funds::parser<hgvsc::reference> ref()
            {
                using namespace varoom::funds;
                return fmap([](accession a) { return hgvsc::reference{a.name, a.version}; }, grammar_basics::accn());
            }
        };

        struct grammar : grammar_basics
        {
            /**
             *
             * exprn <- g_var | c_var | n_var | r_var | p_var
             *
             * g_var <- accn ":g." (g_single | g_multi)
             *
             * <t>_single <- <t>_change
             *
             * <t>_multi <- <t>_phased_alleles | <t>_unphased_alleles
             *
             * <t>_phased_alleles <- <t>_phased_allele (';' <t>_phased_allele)+
             *
             * <t>_unphased_alleles <- <t>_change (';' <t>_change)+
             *
             * <t>_phased_allele <- '[' <t>_single (';' <t>_single)* ']'
             *
             * <t>_change <- <t>_id | <t>_sub | <t>_ins | <t>_del | <t>_delins | <t>_dup | <t>_inv | <t>_con | <t>_rep
             *
             * <t>_id <- <t>_pos '='
             *
             * <t>_sub <- <t>_pos nuc '>' nuc
             *
             * <t>_ins <- <t>_pos '_' <t>_pos "ins" (nuc+|num)
             *
             * <t>_del <- <t>_loc "del" (nuc+|num)?
             *
             * <t>_delins <- <t>_loc "del" (nuc+|num)? "ins" (nuc+|num)
             *
             * <t>_dup <- <t>_loc "dup" (nuc+|num)?
             *
             * <t>_inv <- <t>_loc "inv" (nuc+|num)?
             *
             * <t>_con <- <t>_loc "con" (<t>_loc | <t>_ref_loc)
             *
             * <t>_rep <- <t>_loc (nuc+ '[' num ']')+
             *
             * <t>_ref_loc <- (<t>_ref ':')? <t>_loc
             *
             * <t>_ref <- accn
             *
             * <t>_interval <- <t>_pos '_' <t>_pos?
             *
             * <t>_loc <- <t>_pos ('_' <t>_pos)?      == <t>_interval | <t>_pos
             *
             * g_pos <- num
             *
             * c_pos <- s_pos
             *
             * n_pos <- s_pos
             *
             * s_pos <- ('-' | '*')? num (('+' | '-') num)?
             *
             * refseq <- accn ( '(' accn ')' )?
             *
             */

            static varoom::funds::parser<variant_ptr> g_var()
            {
                using namespace varoom::funds;
                return g_single();
            }

#if 0
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
                return seq<variant_ptr>(sepby1(change(acc), sym(';')), [=](list<variant_ptr> aa) {
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
                return seq<variant_ptr>(sym('['), sepby1(change(acc), sym(';')), sym(']'), [=](char, list<variant_ptr> aa, char) {
                    std::vector<variant_ptr> alleles;
                    while (!aa.empty())
                    {
                        alleles.push_back(aa.head());
                        aa = aa.tail();
                    }
                    return variant_ptr(new hgvsg_allele(acc, alleles, false));
                });
            }
#endif

            static varoom::funds::parser<variant_ptr> g_single()
            {
                using namespace varoom::funds;
                return (seq<genomic_reference>(ref<hgvsg>(), str(":g."),
                    [](genomic_reference acc, std::string) { return acc; }) >>= [](genomic_reference acc) {
                        return change<hgvsg>(acc);
                    });
            }

            template <typename X>
            static varoom::funds::parser<variant_ptr> change(typename X::reference acc)
            {
                using namespace varoom::funds;
                return altx({
                    id<X>(acc), sub<X>(acc),
                    delins<X>(acc),
                    del<X>(acc), ins<X>(acc)
                });
            }

            template <typename X>
            static varoom::funds::parser<variant_ptr> id(typename X::reference acc)
            {
                using locus = typename X::locus;
                using namespace varoom::funds;
                return seq<variant_ptr>(loc<X>(), sym('='), [=](locus loc, char) {
                    return variant_ptr(new hgvs_id<X>(acc, loc));
                });
            }

            template <typename X>
            static varoom::funds::parser<variant_ptr> sub(typename X::reference acc)
            {
                using position = typename X::position;
                using namespace varoom::funds;
                return seq<variant_ptr>(pos<X>(), nuc(), sym('>'), nuc(), [=](position pos, nucleotide ref, char, nucleotide alt) {
                    return variant_ptr(new hgvs_sub<X>(acc, pos, ref, alt));
                });
            }

            template <typename X>
            static varoom::funds::parser<variant_ptr> ins(typename X::reference acc)
            {
                using locus = typename X::locus;
                using namespace varoom::funds;
                return seq<variant_ptr>(interval<X>(), str("ins"), +nuc(), [=](locus loc, std::string, nucleotides alt) {
                    return variant_ptr(new hgvs_ins<X>(acc, loc, alt));
                });
            }

            template <typename X>
            static varoom::funds::parser<variant_ptr> del(typename X::reference acc)
            {
                using locus = typename X::locus;
                using namespace varoom::funds;
                return seq<variant_ptr>(loc<X>(), str("del"), *nuc(), [=](locus loc, std::string, nucleotides ref) {
                    return variant_ptr(new hgvs_del<X>(acc, loc, ref));
                });
            }

            template <typename X>
            static varoom::funds::parser<variant_ptr> delins(typename X::reference acc)
            {
                using locus = typename X::locus;
                using namespace varoom::funds;
                return seq<variant_ptr>(loc<X>(), str("del"), *nuc(), str("ins"), +nuc(),
                    [=](locus loc, std::string, nucleotides ref, std::string, nucleotides alt) {
                        return variant_ptr(new hgvs_delins<X>(acc, loc, ref, alt));
                    });
            }

            template <typename X>
            static varoom::funds::parser<variant_ptr> dup(typename X::reference acc)
            {
                using locus = typename X::locus;
                using namespace varoom::funds;
                return seq<variant_ptr>(loc<X>(), str("dup"), *nuc(), [=](locus loc, std::string, nucleotides ref) {
                    return variant_ptr(new hgvs_dup<X>(acc, loc, ref));
                });
            }

            template <typename X>
            static varoom::funds::parser<variant_ptr> inv(typename X::reference acc)
            {
                using locus = typename X::locus;
                using namespace varoom::funds;
                return seq<variant_ptr>(loc<X>(), str("inv"), *nuc(), [=](locus loc, std::string s, nucleotides ref) {
                    return variant_ptr(new hgvs_inv<X>(acc, loc, ref));
                });
            }

            template <typename X>
            static varoom::funds::parser<variant_ptr> con(typename X::reference acc)
            {
                using locus = typename X::locus;
                using ref_and_loc = typename X::ref_and_loc;
                using namespace varoom::funds;
                return seq<variant_ptr>(loc<X>(), str("con"), ref_loc<X>(acc),
                    [=](locus loc, std::string, ref_and_loc oth) {
                        return variant_ptr(new hgvs_con<X>(acc, loc, oth.first, oth.second));
                });
            }

            template <typename X>
            static varoom::funds::parser<typename X::ref_and_loc> ref_loc(typename X::reference acc)
            {
                using reference = typename X::reference;
                using locus = typename X::locus;
                using ref_and_loc = typename X::ref_and_loc;
                using namespace varoom::funds;

                return altx({
                    seq<ref_and_loc>(ref<X>(), sym(':'), loc<X>(), [](reference ref, char, locus loc) {
                        return ref_and_loc(ref, loc);
                    }),
                    seq<ref_and_loc>(loc<X>(), [=](locus loc) {
                        return ref_and_loc(acc, loc);
                    })
                });
            }

            template <typename X>
            static varoom::funds::parser<typename X::locus> loc()
            {
                using position = typename X::position;
                using locus = typename X::locus;
                using namespace varoom::funds;

                return altx({
                    interval<X>(),
                    fmap([](position p) {
                        return locus(p, p);
                    }, pos<X>())
                });
            }

            template <typename X>
            static varoom::funds::parser<typename X::locus> interval()
            {
                using position = typename X::position;
                using locus = typename X::locus;

                using namespace varoom::funds;
                return seq<locus>(pos<X>(), sym('_'), pos<X>(), [](position f, char c, position l) {
                    return locus(f, l);
                });
            }

            template <typename X>
            static varoom::funds::parser<typename X::position> pos()
            {
                return grammar_traits<X>::pos();
            }

            template <typename X>
            static varoom::funds::parser<typename X::reference> ref()
            {
                return grammar_traits<X>::ref();
            }

        };
    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_HGVS_GRAMMAR_HPP
