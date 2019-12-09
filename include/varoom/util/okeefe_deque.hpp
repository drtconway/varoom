#ifndef VAROOM_UTIL_OKEEFE_DEQUE_HPP
#define VAROOM_UTIL_OKEEFE_DEQUE_HPP

#include <vector>

namespace varoom
{
    template <typename T>
    class okeefe_deque
    {
    public:
        okeefe_deque() {}

        size_t size() const
        {
            return m_front.size() + m_back.size();
        }

        const T& front() const
        {
            if (m_front.size() == 0)
            {
                flip_to_front();
            }
            return m_front.back();
        }

        T& front()
        {
            if (m_front.size() == 0)
            {
                flip_to_front();
            }
            return m_front.back();
        }

        const T& back() const
        {
            if (m_back.size() == 0)
            {
                flip_to_back();
            }
            return m_back.back();
        }

        T& back()
        {
            if (m_back.size() == 0)
            {
                flip_to_back();
            }
            return m_back.back();
        }

        okeefe_deque& push_front(const T& p_item)
        {
            m_front.push_back(p_item);
            return *this;
        }

        okeefe_deque& push_front(T&& p_item)
        {
            m_front.push_back(p_item);
            return *this;
        }

        okeefe_deque& push_back(const T& p_item)
        {
            m_back.push_back(p_item);
            return *this;
        }

        okeefe_deque& push_back(T&& p_item)
        {
            m_back.push_back(p_item);
            return *this;
        }

        okeefe_deque& pop_front()
        {
            if (m_front.size() == 0)
            {
                flip_to_front();
            }
            m_front.pop_back();
            return *this;
        }

        okeefe_deque& pop_back()
        {
            if (m_back.size() == 0)
            {
                flip_to_back();
            }
            m_back.pop_back();
            return *this;
        }

        void swap(okeefe_deque& p_other)
        {
            std::swap(m_front, p_other.m_front);
            std::swap(m_back, p_other.m_back);
        }

    private:
        void flip_to_front()
        {
            while (m_back.size())
            {
                m_front.push_back(std::string());
                std::swap(m_front.back(), m_back.back());
                m_back.pop_back();
            }
        }

        void flip_to_back()
        {
            while (m_front.size())
            {
                m_back.push_back(std::string());
                std::swap(m_front.back(), m_back.back());
                m_front.pop_back();
            }
        }

        std::vector<T> m_front;
        std::vector<T> m_back;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_OKEEFE_DEQUE_HPP
