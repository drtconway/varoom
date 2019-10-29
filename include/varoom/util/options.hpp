#ifndef VAROOM_UTIL_OPTIONS_HPP
#define VAROOM_UTIL_OPTIONS_HPP

#include <functional>
#include <nlohmann/json.hpp>

namespace varoom
{
    namespace util
    {
        typedef std::function<bool(const option_value&)> option_type;
        typedef nlohmann::json option_value;

        class options
        {
        public:
            void required(const std::string& p_name, option_type p_type)
            {
                if (m_required.find(p_name) != m_required.end()
                    || m_optional.find(p_name) != m_optional.end())
                {
                    throw std::logic_error("duplicate required option name");
                }
                m_required[p_name] = p_type;
            }

            void optional(const std::string& p_name, option_type p_type)
            {
                if (m_required.find(p_name) != m_required.end()
                    || m_optional.find(p_name) != m_optional.end())
                {
                    throw std::logic_error("duplicate optional option name");
                }
                m_optional[p_name] = p_type;
            }

            void optional(const std::string& p_name, option_type p_type, const option_value& p_default)
            {
                if (m_required.find(p_name) != m_required.end()
                    || m_optional.find(p_name) != m_optional.end())
                {
                    throw std::logic_error("duplicate optional option name");
                }
                m_optional[p_name] = p_type;
                m_defaults[p_name] = p_default;
            }

        private:
            std::map<std::string,option_type> m_required;
            std::map<std::string,option_type> m_optional;
            std::map<std::string,option_value> m_defaults;
        };

        bool str_option(const option_value& p_value)
        {
            return p_value.is_string();
        }

        bool int_option(const option_value& p_value)
        {
            return p_value.is_number();
        }

        bool flt_option(const option_value& p_value)
        {
            return p_value.is_number();
        }

        bool str_list_option(const option_value& p_value)
        {
            return p_value.is_array();
        }
    }
    // namepsace util
}
// namepsace varoom

#endif // VAROOM_UTIL_OPTIONS_HPP
