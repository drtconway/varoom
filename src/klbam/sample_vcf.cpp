#include "varoom/command.hpp"
#include "varoom/util/typed_tsv.hpp"
#include "varoom/util/files.hpp"
#include "varoom/seq/sequence_factory.hpp"
#include "varoom/vcf/vcf_writer.hpp"
#include "varoom/util/lazy.hpp"

#include <cmath>
#include <initializer_list>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::seq;
using namespace varoom::vcf;

namespace // anonymous
{
    string tabs(initializer_list<const char*> p_parts)
    {
        string s;
        for (auto i = p_parts.begin(); i != p_parts.end(); ++i)
        {
            if (s.size() > 0)
            {
                s.push_back('\t');
            }
            s.insert(s.size(), *i);
        }
        return s;
    }

    typedef pair<string,uint32_t> locus;
    typedef pair<string,size_t> seq_and_count;
    typedef vector<seq_and_count> seq_and_count_list;

    class sample_vcf_command : public varoom::command
    {
    public:
        sample_vcf_command(const string& p_genome_directory,
                            const string& p_input_filename,
                            const string& p_output_filename)
            : m_genome_directory(p_genome_directory),
              m_input_filename(p_input_filename),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            sequence_factory sfac(m_genome_directory);
            vector<string> types{"str", "uint", "uint", "uint", "uint", "uint", "[str=uint]", "flt", "flt"};

            input_file_holder_ptr inp = files::in(m_input_filename);
            for (typed_tsv_reader s(**inp, types); s.more(); ++s)
            {
            }
        }

    private:

        const string m_genome_directory;
        const string m_input_filename;
        const string m_output_filename;
    };

    class sample_vcf_factory : public command_factory
    {
    public:
        sample_vcf_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("convert sample tsv to vcf format");
            opts.add_options()
                ("help,h", "produce help message")
                ("genome,g", po::value<string>()->default_value("."), "genome directory")
                ("input,i", po::value<string>()->default_value("-"), "input filename, defaults to '-' for stdin")
                ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
                ;

            po::positional_options_description pos;
            pos.add("genome", 1);
            pos.add("input", 1);
            pos.add("output", 1);

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

            params["genome"] = vm["genome"].as<string>();
            params["input"] = vm["input"].as<string>();
            params["output"] = vm["output"].as<string>();

            return params;
        }

        virtual command_ptr create(const json& p_params) const
        {
            if (!p_params.is_object())
            {
                return command_ptr();
            }
            string global_fn = p_params["genome"];
            string input_fn = p_params["input"];
            string output_fn = p_params["output"];
            return command_ptr(new sample_vcf_command(global_fn, input_fn, output_fn));
        }
    };
    
    bool reg = command_factory::add("vcf", command_factory_ptr(new sample_vcf_factory));
}
// namespace anonymous


