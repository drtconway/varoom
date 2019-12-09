#ifndef VAROOM_UTIL_LAZY_HPP
#define VAROOM_UTIL_LAZY_HPP

#include <functional>
#include <memory>
#include <mutex>
#include <boost/core/noncopyable.hpp>

namespace varoom
{
    namespace detail
    {
        template <typename T>
        class lazy_thunk : private boost::noncopyable
        {
        public:
            explicit lazy_thunk(std::function<T()> p_func)
                : m_func(p_func)
            {
            }

            T const& get() const
            {
                if (!m_value_ptr.get())
                {
                    std::lock_guard<std::mutex> lk(m_lock);
                    if (!m_value_ptr.get())
                    {
                        m_value_ptr = std::unique_ptr<T>(new T(m_func()));
                    }
                }
                return *m_value_ptr;
            }

        private:
            std::function<T()> m_func;
            mutable std::mutex m_lock;
            mutable std::unique_ptr<T> m_value_ptr;
        };
    }
    // namespace detail

    template <typename T>
    class lazy
    {
    public:
        using thunk = detail::lazy_thunk<T>;

        explicit lazy(std::function<T()> p_func)
            : m_thunk(new thunk(p_func))
        {
        }

        const T& get() const
        {
            return m_thunk->get();
        }

    private:
        std::shared_ptr<thunk> m_thunk;
    };

    template <typename T>
    lazy<T> make_lazy(const T& p_value)
    {
        return lazy<T>([p_value]() { return p_value; });
    }
}
// namespace varoom

#endif // VAROOM_UTIL_LAZY_HPP
