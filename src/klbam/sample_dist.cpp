#include "varoom/command.hpp"
#include "varoom/seq/locus_stream.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/text.hpp"
#include "varoom/util/typed_tsv.hpp"

#include <cmath>
#include <initializer_list>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::seq;

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

    typedef pair<string,size_t> seq_and_count;
    typedef vector<seq_and_count> seq_and_count_list;

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
            output_file_holder_ptr outp = files::out(m_output_filename);
            std::ostream& out = **outp;
            out << text::tabs({"chr", "pos", "ref", "nA", "nC", "nG", "nT", "indel", "kld"}) << std::endl;

            input_file_holder_ptr globp = files::in(m_global_filename);
            input_file_holder_ptr inp = files::in(m_input_filename);

            vector<string> types{"str", "uint", "uint", "uint", "uint", "uint", "[str=uint]"};
            typed_tsv_reader glob(**globp, types);
            typed_tsv_reader sample(**inp, types);

            std::vector<double> gProbs;
            std::vector<double> sProbs;

            locus_id gLoc;
            if (glob.more())
            {
                const typed_tsv_row& gRow = *glob;
                gLoc = locus_id(any_cast<const string&>(gRow[0]), any_cast<uint64_t>(gRow[1]));
            }
            locus_id sLoc;
            if (sample.more())
            {
                const typed_tsv_row& sRow = *sample;
                sLoc = locus_id(any_cast<const string&>(sRow[0]), any_cast<uint64_t>(sRow[1]));
            }

            while (glob.more() && sample.more())
            {
                if (gLoc < sLoc)
                {
                    ++glob;
                    if (glob.more())
                    {
                        const typed_tsv_row& gRow = *glob;
                        gLoc = locus_id(any_cast<const string&>(gRow[0]), any_cast<uint64_t>(gRow[1]));
                    }
                    continue;
                }

                if (sLoc < gLoc)
                {
                    // If the locus isn't in the global table,
                    // then this is the only representative of the locus
                    // so the divergence is zero.
                    //
                    output_row(out, *sample, 0.0);
                    ++sample;
                    if (sample.more())
                    {
                        const typed_tsv_row& sRow = *sample;
                        sLoc = locus_id(any_cast<const string&>(sRow[0]), any_cast<uint64_t>(sRow[1]));
                    }
                    continue;
                }

                const typed_tsv_row& gRow = *glob;
                gProbs.clear();
                gProbs.push_back(any_cast<uint64_t>(gRow[2]));
                gProbs.push_back(any_cast<uint64_t>(gRow[3]));
                gProbs.push_back(any_cast<uint64_t>(gRow[4]));
                gProbs.push_back(any_cast<uint64_t>(gRow[5]));

                const typed_tsv_row& sRow = *sample;
                sProbs.clear();
                sProbs.push_back(any_cast<uint64_t>(sRow[2]));
                sProbs.push_back(any_cast<uint64_t>(sRow[3]));
                sProbs.push_back(any_cast<uint64_t>(sRow[4]));
                sProbs.push_back(any_cast<uint64_t>(sRow[5]));

                const seq_and_count_list& gIndels = any_cast<const seq_and_count_list&>(gRow[6]);
                const seq_and_count_list& sIndels = any_cast<const seq_and_count_list&>(sRow[6]);

                auto gItr = gIndels.begin();
                auto sItr = sIndels.begin();
                while (gItr != gIndels.end() && sItr != sIndels.end())
                {
                    if (*gItr < *sItr)
                    {
                        ++gItr;
                        continue;
                    }

                    if (*sItr < *gItr)
                    {
                        // This shouldn't really happen, in the idealised setting where
                        // all samples are entailed in the global counts, but pragmatically
                        // we can assume that this is the only instance of the indel so the
                        // sample count == global count
                        //
                        gProbs.push_back(sItr->second);
                        sProbs.push_back(sItr->second);
                        ++sItr;
                        continue;
                    }

                    gProbs.push_back(gItr->second);
                    ++gItr;
                    sProbs.push_back(sItr->second);
                    ++sItr;
                }

                // If there are any extra sample counts,
                // they're equivalent to the case above.
                //
                while (sItr != sIndels.end())
                {
                    gProbs.push_back(sItr->second);
                    sProbs.push_back(sItr->second);
                    ++sItr;
                }

                // Finally, to account for the fact that coverage variation
                // does represent real variablility, we include a "null"
                // observation category with a global count of 1 and  a
                // sample count of 1. This vastly reduces the probability
                // of all samples having the same KL-divergence which leads
                // to problems estimating gamma.
                //
                //gProbs.push_back(1);
                //sProbs.push_back(1);

                double gTot = 0;
                double sTot = 0;
                for (size_t i = 0; i < gProbs.size(); ++i)
                {
                    gTot += gProbs[i];
                    sTot += sProbs[i];
                }
                for (size_t i = 0; i < gProbs.size(); ++i)
                {
                    gProbs[i] /= gTot;
                    sProbs[i] /= sTot;
                }

                double kld = kl_divergence(sProbs, gProbs);
                output_row(out, sRow, kld);

                ++glob;
                if (glob.more())
                {
                    const typed_tsv_row& gRow = *glob;
                    gLoc = locus_id(any_cast<const string&>(gRow[0]), any_cast<uint64_t>(gRow[1]));
                }

                ++sample;
                if (sample.more())
                {
                    const typed_tsv_row& sRow = *sample;
                    sLoc = locus_id(any_cast<const string&>(sRow[0]), any_cast<uint64_t>(sRow[1]));
                }
            }

            while (sample.more())
            {
                output_row(out, *sample, 0.0);
                ++sample;
            }
        }

    private:
        static void output_row(ostream& p_out, const typed_tsv_row& p_row, double kld)
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
                  << '\t' << kld
                  << endl;
        }

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

