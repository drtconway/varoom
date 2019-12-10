#ifndef VAROOM_KLBAM_TYPES_HPP
#define VAROOM_KLBAM_TYPES_HPP

#include <initializer_list>
#include <string>
#include <tuple>
#include <boost/flyweight.hpp>

namespace varoom
{
    namespace klbam
    {
        using chrom = boost::flyweight<std::string>;
        using counts_type = std::tuple<chrom,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,std::string>;
        using counts_table = varoom::table::basic_inmemory_table<chrom,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,std::string>;
        using scores_type = std::tuple<double,double>;

        struct counts
        {
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
    }
    // namespace klbam
}
// namespace varoom

#endif // VAROOM_KLBAM_TYPES_HPP
