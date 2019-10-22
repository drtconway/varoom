#ifndef VAROOM_UTIL_SCALAR_HPP
#define VAROOM_UTIL_SCALAR_HPP

#include <cstdint>

namespace varoom
{
    template <typename Tag> class scalar;

    template<typename SrcTag, typename ResTag>
    struct scalar_conversions
    {
        static scalar<ResTag> cast(const scalar<SrcTag>& p_src);
    };

    template <typename Tag>
    class scalar
    {
    public:
        scalar(const std::int64_t& p_value)
            : m_value(p_value)
        {
        }

        std::int64_t operator()() const
        {
            return m_value;
        }

        template <typename ResTag>
        explicit operator scalar<ResTag> () const
        {
            return scalar_conversions<Tag,ResTag>::cast(*this);
        }

    private:
        const std::int64_t m_value;
    };
}
// namespace varoom


#endif // VAROOM_UTIL_SCALAR_HPP
