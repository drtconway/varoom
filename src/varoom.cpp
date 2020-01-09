#include "varoom/command.hpp"
#include "varoom/util/debug.hpp"

#include <iostream>
#include <nlohmann/json.hpp>

using namespace std;
using namespace varoom;
using nlohmann::json;
namespace po = boost::program_options;

namespace // anonymous
{
}
// namespace anonymous

int main(int argc, const char** argv)
{
    po::options_description global("Global options");
    global.add_options()
        ("debug", po::value<string>(), "configure runtime debugging")
        ("help", "show a help message")
        ("command", po::value<string>(), "command to execute")
        ("subargs", po::value<vector<string> >(), "Arguments for command");

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

    if (vm.count("debug"))
    {
        string debug_text = vm["debug"].as<string>();
        json dbgs = json::parse(debug_text);
        if (dbgs.is_string())
        {
            string nm = dbgs;
            if (!debug::exists(nm))
            {
                cerr << "attempt to manipulate non existant debug" << nm << endl;
                return 1;
            }
            debug::get(nm).enable();
        }
        else if (dbgs.is_array())
        {
            for (size_t i = 0; i < dbgs.size(); ++i)
            {
                string nm = dbgs[i];
                if (!debug::exists(nm))
                {
                    cerr << "attempt to manipulate non existant debug" << nm << endl;
                    continue;
                }
                debug::get(nm).enable();
            }
        }
        else if (dbgs.is_object())
        {
            for (json::iterator it = dbgs.begin(); it != dbgs.end(); ++it)
            {
                if (!debug::exists(it.key()))
                {
                    cerr << "attempt to manipulate non existant debug" << it.key() << endl;
                    continue;
                }
                if (it.value().is_boolean())
                {
                    bool v = it.value();
                    debug::get(it.key()).set(v);
                }
                else
                {
                    debug::get(it.key()).data(it.value());
                    debug::get(it.key()).enable();
                }
            }
        }
        else
        {
            cerr << "unexpected JSON: " << dbgs << endl;
            return 1;
        }
    }

    if (vm.count("command"))
    {
        string cmd_name = vm["command"].as<string>();

        if (!command_factory::has_command(cmd_name))
        {
            cerr << global << endl << endl;
            cerr << "unknown command '" << cmd_name << "'." << endl;
            cerr << "to see a list of commands use: varoom --help" << endl;
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
