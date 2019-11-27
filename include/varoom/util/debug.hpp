#ifndef VAROOM_UTIL_DEBUG_HPP
#define VAROOM_UTIL_DEBUG_HPP

#include <map>
#include <vector>
#include <nlohmann/json.hpp>

namespace varoom
{
    class debug
    {
    public:
        debug(const std::string& p_name)
            : m_on(false)
        {
            all()[p_name] = this;
        }

        debug(const std::string& p_name, bool p_on)
            : m_on(p_on)
        {
            all()[p_name] = this;
        }

        bool on() const
        {
            return m_on;
        }

        void enable()
        {
            m_on = true;
        }

        void disable()
        {
            m_on = false;
        }

        void set(bool p_on)
        {
            m_on = p_on;
        }

        const nlohmann::json& data() const
        {
            return m_data;
        }

        void data(const nlohmann::json& p_data)
        {
            m_data = p_data;
        }
            
        static bool exists(const std::string& p_name)
        {
            auto itr = all().find(p_name);
            return itr != all().end();
        }

        static debug& get(const std::string& p_name)
        {
            auto itr = all().find(p_name);
            if (itr == all().end())
            {
                throw std::runtime_error("named debug not found");
            }
            return *(itr->second);
        }

        static std::vector<std::pair<std::string, bool>> list()
        {
            std::vector<std::pair<std::string,bool>> res;
            for (auto i = all().begin(); i != all().end(); ++i)
            {
                res.push_back(std::pair<std::string,bool>(i->first, i->second->on()));
            }
            return res;
        }

    private:
        bool m_on;
        nlohmann::json m_data;

        static std::map<std::string,debug*>& all()
        {
            static std::map<std::string,debug*> s_all;
            return s_all;
        }
    };
}
// namespace varoom

#endif

