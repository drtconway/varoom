#ifndef VAROOM_LAZY_HPP
#define VAROOM_LAZY_HPP

namespace varoom
{
    template <typename T, typename U>
    struct lazy_traits
    {
        static std::shared_ptr<T> make(const U& p_source);
    };

    template <typename T, typename U, typename V = lazy_traits<U, V>>
    class lazy
    {
    public:
        lazy(const U& p_source)
            : m_source(p_source)
        {
        }

        const T& operator*() const
        {
            if (!m_value.get())
            {
                m_value = V::make(m_source);
            }
            return *m_value;
        }

    private:
        const U& m_source;
        mutable std::shared_ptr<T> m_value;
    };
}
// namespace varoom

#endif // VAROOM_LAZY_HPP
