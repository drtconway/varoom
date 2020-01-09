#include "varoom/command.hpp"
#include "varoom/seq/index.hpp"

#include <cmath>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::seq;

namespace // anonymous
{

    using strings = vector<string>;

    class genome_command : public varoom::command
    {
    public:
        genome_command(const strings& p_input_filenames, const string& p_output_filename)
            : m_input_filenames(p_input_filenames),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            varoom::seq::index::make(m_input_filenames, m_output_filename);
        }

    private:
        const strings m_input_filenames;
        const string m_output_filename;
    };

    class genome_factory : public command_factory
    {
    public:
        genome_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("compile a compact genome reference");
            opts.add_options()
                ("help,h", "produce help message")
                ("input,i", po::value<strings>(), "input filename, defaults to '-' for stdin")
                ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
                ;

            po::positional_options_description pos;
            pos.add("output", 1);
            pos.add("input", -1);

            po::variables_map vm;
            po::parsed_options parsed
                = po::command_line_parser(p_args).options(opts).positional(pos).run();
            po::store(parsed, vm);
            po::notify(vm);

            json params;
            if (p_globals.count("help") || vm.count("help"))
            {
                cout << p_global_opts << endl;
                cout << opts << endl;
                return params;
            }

            strings ss;
            if (vm.count("input"))
            {
                ss = vm["input"].as<strings>();
            }
            params["input"] = ss;
            params["output"] = vm["output"].as<string>();

            return params;
        }

        virtual command_ptr create(const json& p_params) const
        {
            if (!p_params.is_object())
            {
                return command_ptr();
            }
            strings input_fn = p_params["input"];
            string output_fn = p_params["output"];
            return command_ptr(new genome_command(input_fn, output_fn));
        }
    };
    
    bool reg = command_factory::add("genome", command_factory_ptr(new genome_factory));
}
// namespace anonymous

