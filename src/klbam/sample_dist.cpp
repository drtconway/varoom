#include "varoom/command.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/text.hpp"
#include "varoom/util/table.hpp"
#include "varoom/klbam/types.hpp"

#include <cmath>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::klbam;

namespace // anonymous
{
    double kl_divergence(const std::vector<double>& P, const std::vector<double>& Q)
    {
        // assert P.size() == Q.size()
        // assert sum(P) == 1.0
        // assert sum(Q) == 1.0
        // assert all(Q_i > 0) unless (P_i == Q_i == 0)

        double d = 0;
        for (size_t i = 0; i < P.size(); ++i)
        {
            if (P[i] == 0)
            {
                continue;
            }
            d += P[i]*std::log(P[i]/Q[i]);
        }
        return d;
    }

    template <int I>
    void add(const counts::tuple& p_lhs, const counts::tuple& p_rhs,
             vector<double>& p_glob, vector<double>& p_sampl, double& p_gt, double& p_st)
    {
        uint32_t g = std::get<I>(p_lhs);
        uint32_t s = std::get<I>(p_rhs);

        p_gt += g;
        p_st += s;

        if (s > 0)
        {
            p_glob.push_back(g);
            p_sampl.push_back(s);
        }
    }

    class sample_dist_command : public varoom::command
    {
    public:
        sample_dist_command(const string& p_global_filename,
                            const string& p_input_filename,
                            const string& p_output_filename)
            : m_global_filename(p_global_filename),
              m_input_filename(p_input_filename),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            map<locus,counts::tuple> glob;
            {
                std::function<void(const counts::tuple&)> f = [&] (const counts::tuple& p_counts) mutable {
                    locus l(std::get<counts::chr>(p_counts), std::get<counts::pos>(p_counts));
                    glob[l] = p_counts;
                };
                input_file_holder_ptr globp = files::in(m_global_filename);
                counts::istream_reader g(**globp, true);
                table_utils::for_each(g, f);
            }

            input_file_holder_ptr inp = files::in(m_input_filename);
            counts::istream_reader in(**inp, true);

            output_file_holder_ptr outp = files::out(m_output_filename);
            scores::ostream_writer out(**outp, scores::labels());

            std::function<void(const counts::tuple&, scores::tuple&)> f = [glob] (const counts::tuple& p_in, scores::tuple& p_out) {

                locus l(std::get<counts::chr>(p_in), std::get<counts::pos>(p_in));

                std::get<scores::chr>(p_out) = l.first;
                std::get<scores::pos>(p_out) = l.second;

                auto gitr = glob.find(l);
                if (gitr == glob.end())
                {
                    std::get<scores::kld>(p_out) = 0.0;
                    std::get<scores::pval>(p_out) = 1.0;
                    return;
                }

                std::vector<double> g;
                double gt = 0;

                std::vector<double> s;
                double st = 0;

                add<counts::nA>(gitr->second, p_in, g, s, gt, st);
                add<counts::nC>(gitr->second, p_in, g, s, gt, st);
                add<counts::nG>(gitr->second, p_in, g, s, gt, st);
                add<counts::nT>(gitr->second, p_in, g, s, gt, st);
                add<counts::nI>(gitr->second, p_in, g, s, gt, st);
                add<counts::nD>(gitr->second, p_in, g, s, gt, st);
                //add<counts::nN>(gitr->second, p_in, g, s);
                //add<counts::nX0>(gitr->second, p_in, p_res);
                //add<counts::nX1>(gitr->second, p_in, p_res);

                for (size_t i = 0; i < g.size(); ++i)
                {
                    g[i] /= gt;
                    s[i] /= st;
                }

                double k = kl_divergence(s, g);

                std::get<scores::kld>(p_out) = k;
                std::get<scores::pval>(p_out) = 1.0;
            };
            table_utils::map(f, in, out);
        }

    private:
        const string m_global_filename;
        const string m_input_filename;
        const string m_output_filename;
    };

    class sample_dist_factory : public command_factory
    {
    public:
        sample_dist_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("compute KL-distances");
            opts.add_options()
                ("help,h", "produce help message")
                ("global,l", po::value<string>(), "global base counts")
                ("input,i", po::value<string>()->default_value("-"), "input filename, defaults to '-' for stdin")
                ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
                ;

            po::positional_options_description pos;
            pos.add("global", 1);
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

            params["global"] = vm["global"].as<string>();
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
            string global_fn = p_params["global"];
            string input_fn = p_params["input"];
            string output_fn = p_params["output"];
            return command_ptr(new sample_dist_command(global_fn, input_fn, output_fn));
        }
    };
    
    bool reg = command_factory::add("dist", command_factory_ptr(new sample_dist_factory));
}
// namespace anonymous

