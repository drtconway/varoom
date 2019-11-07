#include "varoom/command.hpp"
#include "varoom/util/tsv.hpp"
#include "varoom/util/files.hpp"

#include <cmath>
#include <initializer_list>

using namespace std;
using namespace varoom;

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

    struct counts_row
    {
        std::string chr;
        uint32_t    pos;
        std::string ref;
        size_t      nA;
        size_t      nC;
        size_t      nG;
        size_t      nT;

        void fill(const tsv_row& p_row)
        {
            chr = p_row.get<std::string>(0);
            pos = p_row.get<uint32_t>(1);
            ref = p_row.get<std::string>(2);
            nA = p_row.get<size_t>(3);
            nC = p_row.get<size_t>(4);
            nG = p_row.get<size_t>(5);
            nT = p_row.get<size_t>(6);
        }

        bool before(const counts_row& p_other) const
        {
            if (chr == p_other.chr)
            {
                return pos < p_other.pos;
            }
            return chr < p_other.chr;
        }
    };

    double kl_distance(const std::vector<double>& P, const std::vector<double>& Q)
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
            input_file_holder_ptr globp = files::in(m_global_filename);
            input_file_holder_ptr inp = files::in(m_input_filename);
            output_file_holder_ptr outp = files::out(m_output_filename);
            std::ostream& out = **outp;
            out << tabs({"chr", "pos", "ref", "nA", "nC", "nG", "nT", "kld"}) << std::endl;

            tsv_reader glob(**globp);
            tsv_reader sample(**inp);
            counts_row g;
            counts_row s;
            while (glob.more() && sample.more())
            {
                g.fill(*glob);
                s.fill(*sample);
                if (g.before(s))
                {
                    ++glob;
                    continue;
                }
                if (s.before(g))
                {
                    throw std::runtime_error("sample contained position not in global counts");
                }
                const double gTot = g.nA + g.nC + g.nG + g.nT;
                if (gTot == 0)
                {
                    // XXX probably should warn, or something.
                    continue;
                }
                const vector<double> pg{g.nA / gTot, g.nC / gTot, g.nG / gTot, g.nT / gTot};
                const double sTot = s.nA + s.nC + s.nG + s.nT;
                if (sTot == 0)
                {
                    // XXX probably should warn, or something.
                    continue;
                }
                const vector<double> ps{s.nA / sTot, s.nC / sTot, s.nG / sTot, s.nT / sTot};

                double kld = kl_distance(ps, pg);
                out << s.chr 
                    << '\t' << s.pos
                    << '\t' << s.ref
                    << '\t' << s.nA
                    << '\t' << s.nC
                    << '\t' << s.nG
                    << '\t' << s.nT
                    << '\t' << kld
                    << std::endl;
            }
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

            po::variables_map vm;
            po::parsed_options parsed
                = po::command_line_parser(p_args).options(opts).run();
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

