#ifndef VAROOM_UTIL_LAZY_HPP
#define VAROOM_UTIL_LAZY_HPP

#include <functional>

namespace varoom
{
    template <typename T>
    class lazy
    {
    public:
        explicit lazy(std::function<T()> p_func)
            : m_func(p_func), m_thunk(&thunk_force), m_value(T())
        {
        }

        T const& get() const
        {
            return m_thunk(this);
        }
    
    private:
        static T const&  thunk_force(const lazy* p_lazy)
        {
            return p_lazy->set_value();
        }

        static T const& thunk_get(const lazy* p_lazy)
        {
            return p_lazy->get_value();
        }

        T const& get_value() const
        {
            return m_value;
        }

        T const& set_value() const
        {
            m_value = m_func();
            m_thunk = &thunk_get;
            return get_value();
        }

        std::function<T()> m_func;
        mutable T const& (*m_thunk)(const lazy*);
        mutable T m_value;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_LAZY_HPP
