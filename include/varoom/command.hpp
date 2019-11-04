#ifndef VAROOM_COMMAND_HPP
#define VAROOM_COMMAND_HPP

#include <string>
#include <map>
#include <memory>
#include <boost/program_options.hpp>

namespace varoom
{
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

    class command_factory
    {
    public:
        using factory =  command_ptr(*)(const command_parameters& p_params);

        command_factory() = delete;

        static bool has_command(const std::string& p_name)
        {
            auto itr = facs().find(p_name);
            return (itr != facs().end());
        }

        static const command_options& options(const std::string& p_name)
        {
            auto itr = opts().find(p_name);
            if (opts().find(p_name) == opts().end())
            {
                throw std::runtime_error("unknown factory name");
            }
            return *(itr->second);
        }

        static command_ptr create(const std::string& p_name, const command_parameters& p_params)
        {
            auto itr = facs().find(p_name);
            if (itr == facs().end())
            {
                throw std::domain_error("factory name not known");
            }
            return itr->second(p_params);
        }

        static bool add(const std::string& p_name, factory p_factory, command_options_ptr p_options_ptr)
        {
            if (facs().find(p_name) != facs().end())
            {
                // the name is taken!
                return false;
            }
            facs()[p_name] = p_factory;
            opts()[p_name] = p_options_ptr;
            return true;
        }

    private:
        static std::map<std::string,factory>& facs()
        {
            static std::map<std::string,factory> m;
            return m;
        }

        static std::map<std::string,command_options_ptr>& opts()
        {
            static std::map<std::string,command_options_ptr> m;
            return m;
        }
    };


}
// namespace varoom

#endif // VAROOM_COMMAND_HPP
