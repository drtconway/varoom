#include "varoom/command.hpp"

#include <iostream>
#include <unordered_map>
#include <vector>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>
#include <nlohmann/json.hpp>

#include "varoom/kmers.hpp"
#include "varoom/kmers/kmer_accumulator.hpp"
#include "varoom/seq/fastq.hpp"
#include "varoom/seq/fastq_pair.hpp"
#include "varoom/util/files.hpp"

using namespace std;
using namespace boost;
using namespace boost::filesystem;
using namespace varoom;
using namespace varoom::seq;
using nlohmann::json;

namespace // anonymous
{
    typedef vector<string> strings;

    typedef fastq_pair<fastq_reader,fastq_reader> fastq_pair_reader;
    typedef fastq_pair_reader::item_type read_pair;

    class hist_acc
    {
    public:
        hist_acc(unordered_map<size_t,size_t>& p_hist)
            : m_hist(p_hist)
        {
        }

        hist_acc& push_back(const pair<size_t,size_t>& p_item)
        {
            m_hist[p_item.second] += 1;
            return *this;
        }

    private:
        unordered_map<size_t,size_t>& m_hist;
    };

    class hist_command : public varoom::command
    {
    public:
        hist_command(const strings& p_input_filenames,
                            const string& p_output_filename)
            : m_input_filenames(p_input_filenames),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            const size_t K = 21;

            kmer_accumulator X;
            vector<kmer> lhsFwd;
            vector<kmer> lhsRev;
            vector<kmer> rhsFwd;
            vector<kmer> rhsRev;
            for (size_t i = 0; i + 1 < m_input_filenames.size(); i += 2)
            {
                input_file_holder_ptr inp1 = files::in(m_input_filenames[i]);
                input_file_holder_ptr inp2 = files::in(m_input_filenames[i+1]);
                fastq_reader r1(**inp1);
                fastq_reader r2(**inp2);
                for (fastq_pair_reader r1r2(r1, r2); r1r2.more(); ++r1r2)
                {
                    kmers::make(std::get<1>((*r1r2).first), K, lhsFwd, lhsRev);
                    kmers::make(std::get<1>((*r1r2).second), K, rhsFwd, rhsRev);

                    X += lhsFwd;
                    X += lhsRev;
                    X += rhsFwd;
                    X += rhsRev;
                }
            }

            unordered_map<size_t,size_t> H;
            hist_acc h(H);
            X.visit(h);

            for (auto i = H.begin(); i != H.end(); ++i)
            {
                cout << i->first << '\t' << i->second << endl;
            }
        }

    private:

        const strings m_input_filenames;
        const string m_output_filename;
    };

    class hist_factory : public command_factory
    {
    public:
        hist_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const vector<string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("count bases");
            opts.add_options()
                ("help,h", "produce help message")
                ("inputs", po::value<strings>(), "input filenames for paired reads")
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
                params["inputs"] = json::array();
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
            strings inputs_fn = p_params["inputs"];
            string output_fn = p_params["output"];
            return command_ptr(new hist_command(inputs_fn, output_fn));
        }
    };
    
    bool reg = command_factory::add("hist", command_factory_ptr(new hist_factory));
}
// namespace anonymous

