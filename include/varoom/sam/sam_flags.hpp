#ifndef VAROOM_SAM_SAM_FLAGS_HPP
#define VAROOM_SAM_SAM_FLAGS_HPP

namespace varoom
{
    struct sam_flags
    {
        static constexpr uint32_t paired = 1;
        static constexpr uint32_t proper_pair = 2;
        static constexpr uint32_t unmapped = 4;
        static constexpr uint32_t mate_unmapped = 8;
        static constexpr uint32_t reverse = 16;
        static constexpr uint32_t mate_reverse = 32;
        static constexpr uint32_t read_1 = 64;
        static constexpr uint32_t read_2 = 128;
        static constexpr uint32_t secondary = 256;
        static constexpr uint32_t qcfail = 512;
        static constexpr uint32_t dup = 1024;
        static constexpr uint32_t supplementary = 2048;

        static bool is_paired(const uint32_t& p_flgs)
        {
            return p_flgs & paired;
        }

        static bool is_proper_pair(const uint32_t& p_flgs)
        {
            return p_flgs & proper_pair;
        }

        static bool is_unmapped(const uint32_t& p_flgs)
        {
            return p_flgs & unmapped;
        }

        static bool is_mate_unmapped(const uint32_t& p_flgs)
        {
            return p_flgs & mate_unmapped;
        }

        static bool is_reverse(const uint32_t& p_flgs)
        {
            return p_flgs & reverse;
        }

        static bool is_mate_reverse(const uint32_t& p_flgs)
        {
            return p_flgs & mate_reverse;
        }

        static bool is_read_1(const uint32_t& p_flgs)
        {
            return p_flgs & read_1;
        }

        static bool is_read_2(const uint32_t& p_flgs)
        {
            return p_flgs & read_2;
        }

        static bool is_secondary(const uint32_t& p_flgs)
        {
            return p_flgs & secondary;
        }

        static bool is_qcfail(const uint32_t& p_flgs)
        {
            return p_flgs & qcfail;
        }

        static bool is_dup(const uint32_t& p_flgs)
        {
            return p_flgs & dup;
        }

        static bool is_supplementary(const uint32_t& p_flgs)
        {
            return p_flgs & supplementary;
        }
    };
}
// namespace varoom

#endif // VAROOM_SAM_SAM_FLAGS_HPP
