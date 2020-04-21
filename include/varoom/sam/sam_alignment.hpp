#ifndef VAROOM_SAM_SAM_ALIGNMENT_HPP
#define VAROOM_SAM_SAM_ALIGNMENT_HPP

#include <string>

namespace varoom
{
    struct sam_alignment
    {
        std::string name;
        uint32_t flags;
        std::string chr;
        uint32_t pos;
        uint32_t mapq;
        std::string cigar;
        std::string mate_chr;
        uint32_t mate_pos;
        int32_t tlen;
        std::string seq;
        std::string qual;
    };
}
// namespace varoom

#endif // VAROOM_SAM_SAM_ALIGNMENT_HPP
