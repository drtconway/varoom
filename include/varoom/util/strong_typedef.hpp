#ifndef VAROOM_UTIL_STRONG_TYPEDEF_HPP
#define VAROOM_UTIL_STRONG_TYPEDEF_HPP

#include <type_traits>
#include <utility>

namespace varoom
{
    template <typename Tag, typename T>
    class strong_typedef
    {
    public:

        strong_typedef()
            : m_value()
        {
        }

        explicit strong_typedef(const T& p_value)
            : m_value(p_value)
        {
        }

        explicit strong_typedef(T&& p_value) noexcept(std::is_nothrow_move_constructible<T>::value)
            : m_value(std::move(p_value))
        {
        }

        explicit operator T&() noexcept
        {
            return m_value;
        }

        explicit operator const T&() const noexcept
        {
            return m_value;
        }

        friend void swap(strong_typedef& p_lhs, strong_typedef& p_rhs) noexcept
        {
            std::swap(static_cast<T&>(p_lhs), static_cast<T&>(p_rhs));
        }

    private:
        T m_value;
    };

    namespace detail
    {
        template <typename Tag, typename T>
        T underlying_type(strong_typedef<Tag, T>);
    }
    // namespace detail

    template <typename T>
    using underlying_type = decltype(detail::underlying_type(std::declval<typename std::decay<T>::type>()));

    template <class ST>
    struct addition
    {
        friend ST& operator+=(ST& lhs, const ST& rhs)
        {
            using type = underlying_type<ST>;
            static_cast<type&>(lhs) += static_cast<const type&>(rhs);
            return lhs;
        }

        friend ST operator+(const ST& lhs, const ST& rhs)
        {
            using type = underlying_type<ST>;
            return ST(static_cast<const type&>(lhs) + static_cast<const type&>(rhs));
        }
    };

    template <class ST>
    struct subtraction
    {
        friend ST& operator-=(ST& lhs, const ST& rhs)
        {
            using type = underlying_type<ST>;
            static_cast<type&>(lhs) -= static_cast<const type&>(rhs);
            return lhs;
        }

        friend ST operator-(const ST& lhs, const ST& rhs)
        {
            using type = underlying_type<ST>;
            return ST(static_cast<const type&>(lhs) - static_cast<const type&>(rhs));
        }
    };

    template <class ST>
    struct multiplication
    {
        friend ST& operator*=(ST& lhs, const ST& rhs)
        {
            using type = underlying_type<ST>;
            static_cast<type&>(lhs) *= static_cast<const type&>(rhs);
            return lhs;
        }

        friend ST operator*(const ST& lhs, const ST& rhs)
        {
            using type = underlying_type<ST>;
            return ST(static_cast<const type&>(lhs) * static_cast<const type&>(rhs));
        }
    };

    template <class ST>
    struct division
    {
        friend ST& operator/=(ST& lhs, const ST& rhs)
        {
            using type = underlying_type<ST>;
            static_cast<type&>(lhs) /= static_cast<const type&>(rhs);
            return lhs;
        }

        friend ST operator/(const ST& lhs, const ST& rhs)
        {
            using type = underlying_type<ST>;
            return ST(static_cast<const type&>(lhs) / static_cast<const type&>(rhs));
        }
    };

    template <class ST>
    struct modulo
    {
        friend ST& operator%=(ST& lhs, const ST& rhs)
        {
            using type = underlying_type<ST>;
            static_cast<type&>(lhs) %= static_cast<const type&>(rhs);
            return lhs;
        }

        friend ST operator%(const ST& lhs, const ST& rhs)
        {
            using type = underlying_type<ST>;
            return ST(static_cast<const type&>(lhs) % static_cast<const type&>(rhs));
        }
    };

    template <class ST>
    struct integer_arithmetic : addition<ST>, subtraction<ST>, multiplication<ST>, division<ST>, modulo<ST>
    {
    };

}
// namespace varoom

#endif // VAROOM_UTIL_STRONG_TYPEDEF_HPP
