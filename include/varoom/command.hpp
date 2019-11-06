#ifndef VAROOM_COMMAND_HPP
#define VAROOM_COMMAND_HPP

#include <string>
#include <map>
#include <memory>
#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>

namespace varoom
{
    typedef nlohmann::json json;
    typedef boost::program_options::options_description command_options;
    typedef std::shared_ptr<command_options> command_options_ptr;
    typedef boost::program_options::variables_map command_parameters;


    class command
    {
    public:
        virtual void operator()() = 0;

        virtual ~command() {}
    };
    typedef std::shared_ptr<command> command_ptr;

    class command_factory;
    typedef std::shared_ptr<command_factory> command_factory_ptr;

    class command_factory
    {
    public:
        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const = 0;

        virtual command_ptr create(const json& p_params) const = 0;

        static bool has_command(const std::string& p_name)
        {
            auto itr = facs().find(p_name);
            return (itr != facs().end());
        }

        static json parse(const std::string& p_name,
                          const command_options& p_global_opts, const command_parameters& p_globals,
                          const std::vector<std::string>& p_args)
        {
            auto itr = facs().find(p_name);
            if (itr == facs().end())
            {
                throw std::domain_error("factory name not known");
            }
            return itr->second->parse(p_global_opts, p_globals, p_args);
        }

        static command_ptr create(const std::string& p_name, const json& p_params)
        {
            auto itr = facs().find(p_name);
            if (itr == facs().end())
            {
                throw std::domain_error("factory name not known");
            }
            return itr->second->create(p_params);
        }

        static bool add(const std::string& p_name, command_factory_ptr p_factory)
        {
            if (facs().find(p_name) != facs().end())
            {
                // the name is taken!
                return false;
            }
            facs()[p_name] = p_factory;
            return true;
        }

    protected:
        command_factory() {}

    private:
        static std::map<std::string,command_factory_ptr>& facs()
        {
            static std::map<std::string,command_factory_ptr> m;
            return m;
        }
    };


}
// namespace varoom

#endif // VAROOM_COMMAND_HPP
