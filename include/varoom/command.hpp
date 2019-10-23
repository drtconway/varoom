#ifndef VAROOM_COMMAND_HPP
#define VAROOM_COMMAND_HPP

#include <string>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>

namespace varoom
{
    typedef nlohmann::json command_parameters;

    class command
    {
    public:
        virtual void operator()() = 0;

        virtual ~command() {}
    };
    typedef std::shared_ptr<command> command_ptr;

    class command_factory
    {
    public:
        using factory =  command_ptr(*)(const command_parameters& p_params);

        command_factory() = delete;

        static command_ptr create(const std::string& p_name, const command_parameters& p_params)
        {
            auto itr = known().find(p_name);
            if (itr == known().end())
            {
                throw std::domain_error("factory name not known");
            }
            return itr->second(p_params);
        }

        static void add(const std::string& p_name, factory p_factory)
        {
            if (known().find(p_name) != known().end())
            {
                // the name is taken!
                return;
            }
            known()[p_name] = p_factory;
        }

    private:
        static std::map<std::string,factory>& known()
        {
            static std::map<std::string,factory> m;
            return m;
        }
    };


}
// namespace varoom

#endif // VAROOM_COMMAND_HPP
