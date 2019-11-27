#include "varoom/command.hpp"
#include "varoom/seq/locus_stream.hpp"
#include "varoom/util/typed_tsv.hpp"
#include "varoom/util/files.hpp"

#include <cmath>
#include <initializer_list>
#include <boost/math/distributions/gamma.hpp>

using namespace std;
using namespace boost;
using namespace boost::math;
using namespace varoom;
using namespace varoom::seq;

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

    typedef pair<string,size_t> seq_and_count;
    typedef vector<seq_and_count> seq_and_count_list;

    class sample_pval_command : public varoom::command
    {
    public:
        sample_pval_command(const string& p_gamma_filename,
                            const string& p_input_filename,
                            const string& p_output_filename)
            : m_gamma_filename(p_gamma_filename),
              m_input_filename(p_input_filename),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            map<locus_id,pair<double,double>> gammas;
            {
                vector<string> gamma_types{"str", "uint", "flt", "flt"};
                input_file_holder_ptr gammap = files::in(m_gamma_filename);
                for (typed_tsv_reader gam(**gammap, gamma_types); gam.more(); ++gam)
                {
                    const typed_tsv_row& r = *gam;
                    locus_id loc(any_cast<const string&>(r[0]), any_cast<uint64_t>(r[1]));
                    pair<double,double> v(any_cast<double>(r[2]), any_cast<double>(r[3]));
                    gammas[loc] = v;
                }
            }

            output_file_holder_ptr outp = files::out(m_output_filename);
            std::ostream& out = **outp;
            out << tabs({"chr", "pos", "nA", "nC", "nG", "nT", "indel", "kld", "pval"}) << std::endl;

            input_file_holder_ptr inp = files::in(m_input_filename);

            vector<string> types{"str", "uint", "uint", "uint", "uint", "uint", "[str=uint]", "flt"};
            for (typed_tsv_reader sample(**inp, types); sample.more(); ++sample)
            {
                const typed_tsv_row& r = *sample;
                locus_id loc(any_cast<const string&>(r[0]), any_cast<uint64_t>(r[1]));
                const pair<double,double>& gParams = gammas.find(loc)->second;
                gamma_distribution<> G(gParams.first, gParams.second);
                double kld = any_cast<double>(r[7]);
                double pval = cdf(complement(G, kld));
                output_row(out, r, pval);
            }
        }

    private:
        static void output_row(ostream& p_out, const typed_tsv_row& p_row, double p_pval)
        {
            const tsv_column_type& ind = *tsv_column_type::get("[str=uint]");

            string ind_str;
            ind.unmake(p_row[6], ind_str);

            p_out << any_cast<const string&>(p_row[0])
                  << '\t' << any_cast<uint64_t>(p_row[1])
                  << '\t' << any_cast<uint64_t>(p_row[2])
                  << '\t' << any_cast<uint64_t>(p_row[3])
                  << '\t' << any_cast<uint64_t>(p_row[4])
                  << '\t' << any_cast<uint64_t>(p_row[5])
                  << '\t' << ind_str
                  << '\t' << any_cast<double>(p_row[7])
                  << '\t' << p_pval
                  << endl;
        }

        const string m_gamma_filename;
        const string m_input_filename;
        const string m_output_filename;
    };

    class sample_pval_factory : public command_factory
    {
    public:
        sample_pval_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("compute KL-distances");
            opts.add_options()
                ("help,h", "produce help message")
                ("gamma,l", po::value<string>(), "estimated gamma parameters")
                ("input,i", po::value<string>()->default_value("-"), "input filename, defaults to '-' for stdin")
                ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
                ;

            po::positional_options_description pos;
            pos.add("gamma", 1);
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

            params["gamma"] = vm["gamma"].as<string>();
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
            string gamma_fn = p_params["gamma"];
            string input_fn = p_params["input"];
            string output_fn = p_params["output"];
            return command_ptr(new sample_pval_command(gamma_fn, input_fn, output_fn));
        }
    };
    
    bool reg = command_factory::add("pval", command_factory_ptr(new sample_pval_factory));
}
// namespace anonymous

