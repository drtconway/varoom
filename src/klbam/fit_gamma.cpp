#include "varoom/command.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/tsv.hpp"

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

    typedef map<string,size_t> base_counts;
    typedef map<uint32_t,base_counts> pos_base_counts;
    typedef map<string,pos_base_counts> chr_pos_base_counts;

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
            chr_pos_base_counts counts;
            string chr;
            uint32_t pos;
            string seq;
            size_t cnt;
            for (size_t n = 0; n < m_input_filenames.size(); ++n)
            {
                input_file_holder_ptr inp = files::in(m_input_filenames[n]);
                for (auto t = tsv_reader(**inp); t.more(); ++t)
                {
                    const tsv_row& r = *t;
                    chr = r[0];
                    pos = lexical_cast<uint32_t>(make_iterator_range(r[1].first, r[1].second));
                    seq = r[2];
                    cnt = lexical_cast<uint32_t>(make_iterator_range(r[3].first, r[3].second));
                    counts[chr][pos][seq] += cnt;
                }
            }

            output_file_holder_ptr outp = files::out(m_output_filename);
            ostream& out = **outp;
            out << tabs({"chr", "pos", "seq", "count"}) << endl;
            for (auto i = counts.begin(); i != counts.end(); ++i)
            {
                const string& chr = i->first;
                for (auto j = i->second.begin(); j != i->second.end(); ++j)
                {
                    const uint32_t& pos = j->first;
                    for (auto k = j->second.begin(); k != j->second.end(); ++k)
                    {
                        const string& seq = k->first;
                        const size_t& cnt = k->second;
                        out << chr
                            << '\t' << pos
                            << '\t' << seq
                            << '\t' << cnt
                            << endl;
                    }
                }
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

            command_options opts("merge counts");
            opts.add_options()
                ("help,h", "produce help message")
                ("inputs", po::value<strings>(), "input filename, defaults to '-' for stdin")
                ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
                ;

            po::positional_options_description pos;
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
