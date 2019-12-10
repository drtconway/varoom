#ifndef VAROOM_KLBAM_TYPES_HPP
#define VAROOM_KLBAM_TYPES_HPP

#include <initializer_list>
#include <string>
#include <tuple>
#include <boost/flyweight.hpp>

#ifndef VAROOM_SEQ_LOCUS_ORDERING_HPP
#include "varoom/seq/locus_ordering.hpp"
#endif

namespace varoom
{
    namespace klbam
    {
        struct chrom : boost::flyweight<std::string>
        {
            chrom() {}

            chrom(const std::string& p_str)
                : boost::flyweight<std::string>(p_str)
            {
            }

            bool operator<(const chrom& p_other) const
            {
                return varoom::locus_ordering::less(*this, p_other);
            }
        };

        struct counts :  varoom::table<chrom,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,std::string>
        {
            using tuple = tuple_type;
            using table = basic_inmemory_table;
            using istream_reader = basic_istream_table;
            using ostream_writer = basic_ostream_table;
            using read_iterator = basic_read_iterator;
            using write_iterator = basic_write_iterator;

            static constexpr std::size_t chr = 0;
            static constexpr std::size_t pos = 1;
            static constexpr std::size_t nA  = 2;
            static constexpr std::size_t nC  = 3;
            static constexpr std::size_t nG  = 4;
            static constexpr std::size_t nT  = 5;
            static constexpr std::size_t nN  = 6;
            static constexpr std::size_t nX0 = 7;
            static constexpr std::size_t nX1 = 8;
            static constexpr std::size_t nI  = 9;
            static constexpr std::size_t nD  = 10;
            static constexpr std::size_t oth = 11;

            static std::initializer_list<std::string> labels()
            {
                return {"chr", "pos", "nA", "nC", "nG", "nT", "nN", "nX0", "nX1", "nI", "nD", "indels"};
            }
        };

        using scores_table_type = varoom::table<chrom,uint32_t,double,double>;

    }
    // namespace klbam
}
// namespace varoom

#endif // VAROOM_KLBAM_TYPES_HPP
