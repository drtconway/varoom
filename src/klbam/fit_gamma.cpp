#include "varoom/command.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/typed_tsv.hpp"

#include <cmath>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;
using namespace varoom;

namespace // anonymous
{
    typedef vector<string> strings;

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

    class fit_gamma_command : public varoom::command
    {
    public:
        fit_gamma_command(const strings& p_input_filenames,
                            const string& p_output_filename)
            : m_input_filenames(p_input_filenames),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            vector<string> types{"str", "uint", "uint", "uint", "uint", "uint", "[str=uint]", "flt"};

            map<locus,vector<double>> klds;
            for (size_t n = 0; n < m_input_filenames.size(); ++n)
            {
                input_file_holder_ptr inp = files::in(m_input_filenames[n]);
                for (auto t = typed_tsv_reader(**inp, types); t.more(); ++t)
                {
                    const typed_tsv_row& r = *t;
                    locus loc(any_cast<const string&>(r[0]), any_cast<uint64_t>(r[1]));
                    double kld = any_cast<double>(r[7]);
                    klds[loc].push_back(kld);
                }
            }

            output_file_holder_ptr outp = files::out(m_output_filename);
            ostream& out = **outp;
            out << tabs({"chr", "pos", "k", "theta"}) << endl;
            for (auto itr = klds.begin(); itr != klds.end(); ++itr)
            {
                const locus& loc = itr->first;
                const vector<double>& xs = itr->second;
                double n = xs.size();
                double sx = 0.0;
                double sxl = 0.0;
                double sl = 0.0;
                for (size_t i = 0; i < xs.size(); ++i)
                {
                    double x = xs[i];
                    double lx = std::log(x);
                    sx += x;
                    sl += lx;
                    sxl += x*lx;
                }
                double k = n*sx / (n*sxl - sl*sx);
                double theta = (n*sxl - sl*sx)/(n*n);
                out << loc.first
                    << '\t' << loc.second
                    << '\t' << k
                    << '\t' << theta
                    << endl;
            }
        }

    private:
        const strings m_input_filenames;
        const string m_output_filename;
    };

    class fit_gamma_factory : public command_factory
    {
    public:
        fit_gamma_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("fit per locus distributions");
            opts.add_options()
                ("help,h", "produce help message")
                ("inputs", po::value<strings>(), "input filename, defaults to '-' for stdin")
                ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
                ;

            po::positional_options_description pos;
            pos.add("inputs", -1);

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

            if (vm.count("inputs"))
            {
                params["inputs"] = vm["inputs"].as<strings>();
            }
            else
            {
                strings ss{"-"};
                params["inputs"] = ss;
            }
            params["output"] = vm["output"].as<string>();

            return params;
        }

        virtual command_ptr create(const json& p_params) const
        {
            if (!p_params.is_object())
            {
                return command_ptr();
            }
            strings input_fns = p_params["inputs"];
            string output_fn = p_params["output"];
            return command_ptr(new fit_gamma_command(input_fns, output_fn));
        }
    };

    bool reg = command_factory::add("fit", command_factory_ptr(new fit_gamma_factory));
}
// namespace anonymous
