#ifndef VAROOM_HGVS_VARIANT_HPP
#define VAROOM_HGVS_VARIANT_HPP

#include <memory>
#include <string>
#include <vector>

#ifndef VAROOM_UTIL_STRONG_TYPEDEF_HPP
#include "varoom/util/strong_typedef.hpp"
#endif

namespace varoom
{
    namespace hgvs
    {
        struct position
        {
            virtual ~position() {}
        };

        struct genomic_position : position
        {
        };
        using genomic_position_ptr = std::shared_ptr<genomic_position>;

        struct known_genomic_position : genomic_position, strong_typedef<known_genomic_position, uint64_t>
        {
            using strong_typedef::strong_typedef;
        };

        struct unknown_genomic_position : genomic_position { };
        
        struct uncertain_genomic_position : genomic_position
        {
            genomic_position_ptr first;
            genomic_position_ptr last;

            uncertain_genomic_position(genomic_position_ptr fst, genomic_position_ptr lst)
                : first(fst), last(lst)
            {
            }
        };

        struct genomic_reference
        {
            std::string name;
            uint64_t version;
        };

        enum tx_section { UTR5, TR, UTR3 };
        struct spliced_position
        {
            tx_section sec;
            int64_t pos;
            int64_t off;
        };

        using nucleotide = char;
        using nucleotides = std::string;

        struct variant
        {
            virtual ~variant() {}
        };
        using variant_ptr = std::shared_ptr<variant>;

        struct hgvsg : variant
        {
            using reference = genomic_reference;
            using position = uint64_t;
            using locus = std::pair<position,position>;
            using ref_and_loc = std::pair<reference,locus>;
            static constexpr char marker = 'g';

            using general_position = genomic_position_ptr;

            static general_position make_known(position p)
            {
                return general_position(new known_genomic_position(p));
            }

            static general_position make_unknown()
            {
                return general_position(new unknown_genomic_position);
            }

            static general_position make_uncertain(general_position fst, general_position lst)
            {
                return general_position(new uncertain_genomic_position(fst, lst));
            }
        };

        struct hgvsc : variant
        {
            using reference = genomic_reference;
            using position = spliced_position;
            using locus = std::pair<position,position>;
            using ref_and_loc = std::pair<reference,locus>;
            static constexpr char marker = 'c';
        };

        template <typename X>
        struct hgvs_id : X
        {
            using reference = typename X::reference;
            using locus = typename X::locus;

            hgvs_id(const reference& p_acc, const locus& p_loc)
                : acc(p_acc), loc(p_loc)
            {
            }

            reference acc;
            locus loc;
        };

        template <typename X>
        struct hgvs_sub : X
        {
            using reference = typename X::reference;
            using position = typename X::position;

            hgvs_sub(const reference& p_acc, const position& p_pos, const nucleotide& p_ref, const nucleotide& p_alt)
                : acc(p_acc), pos(p_pos), ref(p_ref), alt(p_alt)
            {
            }

            reference acc;
            position pos;
            nucleotide ref;
            nucleotide alt;
        };

        template <typename X>
        struct hgvs_ins : X
        {
            using reference = typename X::reference;
            using locus = typename X::locus;

            hgvs_ins(const reference& p_acc, const locus& p_loc, const nucleotides& p_alt)
                : acc(p_acc), loc(p_loc), alt(p_alt)
            {
            }

            reference acc;
            locus loc;
            nucleotides alt;
        };

        template <typename X>
        struct hgvs_del : X
        {
            using reference = typename X::reference;
            using locus = typename X::locus;

            hgvs_del(const reference& p_acc, const locus& p_loc, const nucleotides& p_ref)
                : acc(p_acc), loc(p_loc), ref(p_ref)
            {
            }

            reference acc;
            locus loc;
            nucleotides ref;
        };

        template <typename X>
        struct hgvs_delins : X
        {
            using reference = typename X::reference;
            using locus = typename X::locus;

            hgvs_delins(const reference& p_acc, const locus& p_loc, const nucleotides& p_ref, const nucleotides& p_alt)
                : acc(p_acc), loc(p_loc), ref(p_ref), alt(p_alt)
            {
            }

            reference acc;
            locus loc;
            nucleotides ref;
            nucleotides alt;
        };

        template <typename X>
        struct hgvs_dup : X
        {
            using reference = typename X::reference;
            using locus = typename X::locus;

            hgvs_dup(const reference& p_acc, const locus& p_loc, const nucleotides& p_ref)
                : acc(p_acc), loc(p_loc), ref(p_ref)
            {
            }

            reference acc;
            locus loc;
            nucleotides ref;
        };

        template <typename X>
        struct hgvs_inv : X
        {
            using reference = typename X::reference;
            using locus = typename X::locus;

            hgvs_inv(const reference& p_acc, const locus& p_loc, const nucleotides& p_ref)
                : acc(p_acc), loc(p_loc), ref(p_ref)
            {
            }

            reference acc;
            locus loc;
            nucleotides ref;
        };

        template <typename X>
        struct hgvs_con : X
        {
            using reference = typename X::reference;
            using locus = typename X::locus;

            hgvs_con(const reference& p_acc, const locus& p_loc, const reference& p_other_acc, const locus& p_other_loc)
                : acc(p_acc), loc(p_loc), other_acc(p_other_acc), other_loc(p_other_loc)
            {
            }

            reference acc;
            locus loc;
            reference other_acc;
            locus other_loc;
        };

    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_VARIANT_HPP
