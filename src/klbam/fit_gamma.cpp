#include "varoom/command.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/table.hpp"
#include "varoom/util/gamma_estimator.hpp"
#include "varoom/klbam/types.hpp"

#include <cmath>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::klbam;

namespace // anonymous
{
    typedef vector<string> strings;


    class fit_gamma_command : public varoom::command
    {
    public:
        fit_gamma_command(bool p_do_merge,
                          const strings& p_input_filenames,
                          const string& p_output_filename)
            : m_do_merge(p_do_merge),
              m_input_filenames(p_input_filenames),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            map<locus,gamma_estimator> klds;
            if (!m_do_merge)
            {
                // Process KL-Divergence scores to produce estimates
                //
                std::function<void(const scores::tuple&)> f = [&] (const scores::tuple& p_item) mutable {
                    locus l(std::get<scores::chr>(p_item), std::get<scores::pos>(p_item));
                    klds[l].push_back(std::get<scores::kld>(p_item));
                };

                for (size_t n = 0; n < m_input_filenames.size(); ++n)
                {
                    input_file_holder_ptr inp = files::in(m_input_filenames[n]);
                    scores::istream_reader in(**inp);
                    table_utils::for_each(in, f);
                }
            }
            else
            {
                // Merge estimates
                //
                std::function<void(const gamma::tuple&)> f = [&] (const gamma::tuple& p_item) mutable {
                    locus l(std::get<gamma::chr>(p_item), std::get<gamma::pos>(p_item));
                    klds[l] += gamma_estimator(std::get<gamma::n>(p_item),
                                               std::get<gamma::sx>(p_item), std::get<gamma::slx>(p_item),
                                               std::get<gamma::sxlx>(p_item));
                };

                for (size_t n = 0; n < m_input_filenames.size(); ++n)
                {
                    input_file_holder_ptr inp = files::in(m_input_filenames[n]);
                    gamma::istream_reader in(**inp);
                    table_utils::for_each(in, f);
                }
            }
            output_file_holder_ptr outp = files::out(m_output_filename);
            gamma::ostream_writer out(**outp, gamma::labels());
            for (auto itr = klds.begin(); itr != klds.end(); ++itr)
            {
                const locus& l = itr->first;
                const gamma_estimator& g = itr->second;
                pair<double,double> v = g();

                gamma::tuple t(l.first, l.second, g.count(), g.sum_x(), g.sum_log_x(), g.sum_x_log_x(), v.first, v.second);
                out << t;
            }
        }

    private:
        const bool m_do_merge;
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
                ("merge,m", "merge estimates, rather than constructing estimates from divergences")
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

            bool mrg = (vm.count("merge") != 0);
            params["merge"] = mrg;

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
            bool do_merge = p_params["merge"];
            strings input_fns = p_params["inputs"];
            string output_fn = p_params["output"];
            return command_ptr(new fit_gamma_command(do_merge, input_fns, output_fn));
        }
    };

    bool reg = command_factory::add("fit", command_factory_ptr(new fit_gamma_factory));
}
// namespace anonymous
