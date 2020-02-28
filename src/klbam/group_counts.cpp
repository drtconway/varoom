#include "varoom/command.hpp"
#include "varoom/util/ranges.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/tsv.hpp"
#include "varoom/util/table.hpp"
#include "varoom/klbam/types.hpp"

#include <unordered_set>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::klbam;

namespace // anonymous
{
    using gene_name = flyweight<string>;

    class annotation_index
    {
    public:
        annotation_index(const string& p_filename)
        {
            map<string,ranges_builder> bldrs;

            input_file_holder_ptr inp = files::in(p_filename);
            for (tsv_reader in(**inp); in.more(); ++in)
            {
                const tsv_row& r = *in;

                string chr = string(r[0]);
                ranges_builder::range v(r.get<size_t>(1), r.get<size_t>(2));
                bldrs[chr].insert(v);
                m_index[chr][v] = string(r[3]);
            }

            for (auto itr = bldrs.begin(); itr != bldrs.end(); ++itr)
            {
                m_ranges[itr->first].make(itr->second);
                //std::cerr << itr->first << '\t' << itr->second.size() << std::endl;
            }
        }

        void operator()(const string& p_chr, const size_t& p_pos, unordered_set<gene_name>& p_res) const
        {
            if (m_ranges.find(p_chr) == m_ranges.end())
            {
                return;
            }
            vector<ranges::range> r;
            m_ranges.at(p_chr).ranges_at(p_pos, r);
            for (auto itr = r.begin(); itr != r.end(); ++itr)
            {
                p_res.insert(m_index.at(p_chr).at(*itr));
            }
            if (0 && p_res.size())
            {
                std::cerr << p_chr << '\t' << p_pos;
                for (auto itr = p_res.begin(); itr != p_res.end(); ++itr)
                {
                    std::cerr << '\t' << (*itr);
                }
                std::cerr << std::endl;
            }
        }

    private:
        map<string,ranges> m_ranges;
        map<string,map<ranges::range,gene_name>> m_index;
    };

    class group_counts_command : public varoom::command
    {
    public:
        group_counts_command(const string& p_regions_filename,
                            const string& p_input_filename,
                            const string& p_output_filename)
            : m_input_filename(p_input_filename),
              m_output_filename(p_output_filename),
              m_regions_index(p_regions_filename)
        {
        }

        virtual void operator()()
        {
            map<gene_name,pair<size_t,size_t>> res;
            {
                std::function<void(const counts::tuple&)> f = [&] (const counts::tuple& p_counts) mutable {
                    string chr = std::get<counts::chr>(p_counts);
                    size_t pos = std::get<counts::pos>(p_counts);
                    size_t cnt =  0;
                    cnt += std::get<counts::nA>(p_counts);
                    cnt += std::get<counts::nC>(p_counts);
                    cnt += std::get<counts::nG>(p_counts);
                    cnt += std::get<counts::nT>(p_counts);
                    cnt += std::get<counts::nI>(p_counts);
                    cnt += std::get<counts::nD>(p_counts);
                    unordered_set<gene_name> gs;
                    m_regions_index(chr, pos, gs);
                    for (auto itr = gs.begin(); itr != gs.end(); ++itr)
                    {
                        res[*itr].first += 1;
                        res[*itr].second += cnt;
                    }
                };
                input_file_holder_ptr inp = files::in(m_input_filename);
                counts::istream_reader in(**inp, true);
                table_utils::for_each(in, f);

                //std::cerr << "finished reading" << std::endl;
            }

            //std::cerr << "opening " << m_output_filename << std::endl;
            output_file_holder_ptr outp = files::out(m_output_filename);
            genes::ostream_writer out(**outp, genes::labels());
            //std::cerr << "writing " << res.size() << " items." << std::endl;
            for (auto itr = res.begin(); itr != res.end(); ++itr)
            {
                //std::cerr << itr->first << '\t' << itr->second.first << '\t' << itr->second.second << std::endl;
                genes::tuple r;
                std::get<genes::gene>(r) = itr->first;
                std::get<genes::len>(r) = itr->second.first;
                std::get<genes::cnt>(r) = itr->second.second;
                out << r;
            }
        }

    private:
        const string m_input_filename;
        const string m_output_filename;
        const annotation_index m_regions_index;
    };

    class group_counts_factory : public command_factory
    {
    public:
        group_counts_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const vector<string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("group counts according to annotations");
            opts.add_options()
                ("help,h", "produce help message")
                ("regions,r", po::value<string>(), "only include counts from regions in the named BED file")
                ("input,i", po::value<string>()->default_value("-"), "input filename, defaults to '-' for stdin")
                ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
                ;

            po::positional_options_description pos;
            pos.add("regions", 1);
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

            params["regions"] = vm["regions"].as<string>();
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
            string regions_fn = p_params["regions"];
            string input_fn = p_params["input"];
            string output_fn = p_params["output"];
            return command_ptr(new group_counts_command(regions_fn, input_fn, output_fn));
        }
    };
    
    bool reg = command_factory::add("group", command_factory_ptr(new group_counts_factory));
}
// namespace anonymous

