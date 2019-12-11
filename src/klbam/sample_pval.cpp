#include "varoom/command.hpp"
#include "varoom/seq/locus_stream.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/gamma_estimator.hpp"
#include "varoom/util/table.hpp"
#include "varoom/klbam/types.hpp"

#include <cmath>
#include <initializer_list>
#include <boost/math/distributions/gamma.hpp>

using namespace std;
using namespace boost;
using namespace boost::math;
using namespace varoom;
using namespace varoom::klbam;

namespace // anonymous
{
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
            map<locus,pair<double,double>> gammas;
            {
                std::function<void(const gamma::tuple&)> f = [&] (const gamma::tuple& p_item) mutable {
                    locus l(std::get<gamma::chr>(p_item), std::get<gamma::pos>(p_item));
                    gamma_estimator g(std::get<gamma::n>(p_item),
                                      std::get<gamma::sx>(p_item), std::get<gamma::slx>(p_item),
                                      std::get<gamma::sxlx>(p_item));
                    gammas[l] = g();
                };
                input_file_holder_ptr inp = files::in(m_gamma_filename);
                gamma::istream_reader in(**inp, true);
                table_utils::for_each(in, f);
            }

            std::function<void(const scores::tuple&, scores::tuple&)> f = [&] (const scores::tuple& p_x, scores::tuple& p_y) {
                locus l(std::get<scores::chr>(p_x), std::get<scores::pos>(p_x));

                auto g_itr = gammas.find(l);
                if (g_itr == gammas.end())
                {
                    throw std::runtime_error("locus not found");
                }
                const pair<double,double>& g = g_itr->second;
                gamma_distribution<> G(g.first, g.second);

                double kld = std::get<scores::kld>(p_x);
                double pval = cdf(complement(G, kld));

                std::get<scores::chr>(p_y) = l.first;
                std::get<scores::pos>(p_y) = l.second;
                std::get<scores::kld>(p_y) = kld;
                std::get<scores::pval>(p_y) = pval;
            };

            output_file_holder_ptr outp = files::out(m_output_filename);
            scores::ostream_writer out(**outp, scores::labels());

            input_file_holder_ptr inp = files::in(m_input_filename);
            scores::istream_reader in(**inp, true);
            table_utils::map(f, in, out);
        }

    private:
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

