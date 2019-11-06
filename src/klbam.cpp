#include "varoom/command.hpp"

#include <iostream>

using namespace std;
using namespace varoom;
namespace po = boost::program_options;

namespace // anonymous
{

}
// namespace anonymous

int main(int argc, const char** argv)
{
    po::options_description global("Global options");
    global.add_options()
        ("debug", "Turn on debug output")
        ("help", "show a help message")
        ("command", po::value<std::string>(), "command to execute")
        ("subargs", po::value<std::vector<std::string> >(), "Arguments for command");

    po::positional_options_description pos;
    pos.add("command", 1).
        add("subargs", -1);

    po::variables_map vm;

    po::parsed_options parsed = po::command_line_parser(argc, argv).
        options(global).
        positional(pos).
        allow_unregistered().
        run();

    po::store(parsed, vm);
    po::notify(vm);

    if (vm.count("command"))
    {
        std::string cmd_name = vm["command"].as<std::string>();

        if (!command_factory::has_command(cmd_name))
        {
            cerr << global << endl << endl;
            cerr << "unknown command '" << cmd_name << "'." << endl;
            cerr << "to see a list of commands use: klbam --help" << endl;
            return 1;
        }

        vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());

        json args = command_factory::parse(cmd_name, global, vm, opts);
        command_ptr cmdp = command_factory::create(cmd_name, args);
        if (cmdp)
        {
            (*cmdp)();
        }
    }

    if (vm.count("help"))
    {
        cout << global << endl;
        return 0;
    }

    return 0;
}
