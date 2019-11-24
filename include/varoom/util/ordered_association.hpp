#ifndef VAROOM_UTIL_ORDERED_ASSOCIATION_HPP
#define VAROOM_UTIL_ORDERED_ASSOCIATION_HPP

#include <functional>
#include <utility>

namespace varoom
{
    namespace detail
    {
        std::string& a_string_reference()
        {
            static std::string s;
            return s;
        }
    }
    // namespace detail

    class ordered_association
    {
    public:
        using const_pair_type = std::pair<const std::string&,const std::string&>;

        ordered_association()
        {
        }

        size_t size() const
        {
            return m_keys.size();
        }

        const_pair_type operator[](size_t p_idx) const
        {
            return const_pair_type(m_keys[p_idx], m_values[p_idx]);
        }

        const std::string& operator[](const std::string& p_key) const
        {
            int i = find(p_key);
            if (i < 0)
            {
                throw std::runtime_error("key not found");
            }
            return m_values[i];
        }

        std::string& operator[](const std::string& p_key)
        {
            int i = find(p_key);
            if (i < 0)
            {
                i = m_keys.size();
                m_keys.push_back(p_key);
                m_values.push_back(std::string());
            }
            return m_values[i];
        }

        void clear()
        {
            m_keys.clear();
            m_values.clear();
        }

    private:

        int find(const std::string& p_key) const
        {
            for (size_t i = 0; i < m_keys.size(); ++i)
            {
                if (m_keys[i] == p_key)
                {
                    return i;
                }
            }
            return -1;
        }

        std::vector<std::string> m_keys;
        std::vector<std::string> m_values;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_ORDERED_ASSOCIATION_HPP
