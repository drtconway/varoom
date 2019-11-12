#include "varoom/command.hpp"
#include "varoom/util/typed_tsv.hpp"
#include "varoom/util/files.hpp"
#include "varoom/seq/sequence_factory.hpp"
#include "varoom/vcf/vcf_writer.hpp"
#include "varoom/util/lazy.hpp"
#include "varoom/util/text.hpp"

#include <cmath>
#include <initializer_list>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::seq;
using namespace varoom::vcf;

namespace // anonymous
{
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
            static uint64_t minCov = 100;
            static uint64_t minAlleleCov = 3;
            static double pval_cutoff = 1e-4;

            sequence_factory sfac(m_genome_directory);

            output_file_holder_ptr outp = files::out(m_output_filename);
            vcf_writer out(**outp);

            vector<string> types{"str", "uint", "uint", "uint", "uint", "uint", "[str=uint]", "flt", "flt"};
            input_file_holder_ptr inp = files::in(m_input_filename);
            for (typed_tsv_reader s(**inp, types); s.more(); ++s)
            {
                const typed_tsv_row& r = *s;
                locus loc(any_cast<const string&>(r[0]), any_cast<uint64_t>(r[1]));
                uint64_t nA = any_cast<uint64_t>(r[2]);
                uint64_t nC = any_cast<uint64_t>(r[3]);
                uint64_t nG = any_cast<uint64_t>(r[4]);
                uint64_t nT = any_cast<uint64_t>(r[5]);
                uint64_t cov = nA + nC + nG + nT;
                const seq_and_count_list& indels = any_cast<const seq_and_count_list&>(r[6]);
                double kld = any_cast<double>(r[7]);
                double pval = any_cast<double>(r[8]);

                if (cov < minCov || pval > pval_cutoff)
                {
                    continue;
                }

                const std::string& seq = sfac[loc.first].second;
                const char ref = text::to_upper(seq[loc.second - 1]); // Strings are 0-indexed.

                if (nA >= minAlleleCov && ref != 'A')
                {
                    make_vcf_row(out, loc, ref, 'A', nA, cov, kld, pval);
                }
                if (nC >= minAlleleCov && ref != 'C')
                {
                    make_vcf_row(out, loc, ref, 'C', nC, cov, kld, pval);
                }
                if (nG >= minAlleleCov && ref != 'G')
                {
                    make_vcf_row(out, loc, ref, 'G', nG, cov, kld, pval);
                }
                if (nT >= minAlleleCov && ref != 'T')
                {
                    make_vcf_row(out, loc, ref, 'T', nT, cov, kld, pval);
                }
                // XXX indels
            }
        }

    private:
        void make_vcf_row(vcf_writer& p_out, const locus& p_loc, char p_ref, char p_alt, uint64_t p_altCov, uint64_t p_cov, double p_kld, double p_pval)
        {
            static std::string dot(".");
            static std::string pass("PASS");
            double phred = -10*std::log10(1 - p_pval);
            vcf_info ifo;
            p_out(p_loc.first, p_loc.second, dot, std::string(1, p_ref),  std::string(1, p_alt), phred, pass, make_lazy(ifo), make_lazy(std::vector<vcf_info>()));
        }

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


