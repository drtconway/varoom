#include "varoom/command.hpp"

#include "varoom/seq/fastq.hpp"
#include "varoom/util/files.hpp"

#include <boost/format.hpp>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::seq;

namespace // anonymous
{
    using strings = vector<string>;

    class scatter_command : public varoom::command
    {
    public:
        scatter_command(const string& p_input_filename, const strings& p_output_filenames)
            : m_input_filename(p_input_filename),
              m_output_filenames(p_output_filenames)
        {
        }

        virtual void operator()()
        {
            vector<output_file_holder_ptr> outs;
            for (size_t i = 0; i < m_output_filenames.size(); ++i)
            {
                outs.push_back(files::out(m_output_filenames[i]));
            }
            const size_t N = outs.size();

            input_file_holder_ptr inp = files::in(m_input_filename);

            size_t n = 0;
            for (fastq_reader r(**inp); r.more(); ++r, ++n)
            {
                while (n >= N)
                {
                    n -= N;
                }
                fastq_writer::write(**(outs[n]), *r);
            }
        }

    private:
        const string m_input_filename;
        const strings m_output_filenames;
    };

    class scatter_factory : public command_factory
    {
    public:
        scatter_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("compile a compact genome reference");
            opts.add_options()
                ("help,h", "produce help message")
                ("input,i", po::value<string>()->default_value("-"), "input filename, defaults to '-' for stdin")
                ("directory,d", po::value<string>(), "directory for output files")
                ("number,n", po::value<int>()->default_value(2), "number of files to split the output into")
                ;

            po::positional_options_description pos;

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

            int n = vm["number"].as<int>();
            string d = vm["directory"].as<string>();

            strings outputs;
            for (int i = 0; i < n; ++i)
            {
                outputs.push_back(str(format("%s/part-%03d.fastq.gz") % d % i));
            }

            params["input"] = vm["input"].as<string>();
            //params["output"] = vm["output"].as<strings>();
            params["output"] = outputs;

            return params;
        }

        virtual command_ptr create(const json& p_params) const
        {
            if (!p_params.is_object())
            {
                return command_ptr();
            }
            string input_fn = p_params["input"];
            strings output_fns = p_params["output"];
            return command_ptr(new scatter_command(input_fn, output_fns));
        }
    };
    
    bool reg = command_factory::add("scatter", command_factory_ptr(new scatter_factory));
}
// namespace anonymous

