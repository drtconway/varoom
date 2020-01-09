#include "varoom/command.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/ranges.hpp"
#include "varoom/util/table.hpp"

#include <cmath>
#include <map>
#include <unordered_map>
#include <boost/flyweight.hpp>

using namespace std;
using namespace boost;
using namespace varoom;

namespace // anonymous
{

    using strings = vector<string>;

    using atom = boost::flyweight<std::string>;

    struct refgene :  varoom::table<uint32_t,std::string,atom,atom,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,std::string,std::string,uint32_t,std::string,atom,atom,std::string>
    {
        using tuple = tuple_type;
        using table = basic_inmemory_table;
        using istream_reader = basic_istream_table;
        using ostream_writer = basic_ostream_table;
        using read_iterator = basic_read_iterator;
        using write_iterator = basic_write_iterator;

        static constexpr std::size_t bin   = 0;
        static constexpr std::size_t acc   = 1;
        static constexpr std::size_t chr   = 2;
        static constexpr std::size_t str   = 3;
        static constexpr std::size_t txs   = 4;
        static constexpr std::size_t txe   = 5;
        static constexpr std::size_t cdss  = 6;
        static constexpr std::size_t cdse  = 7;
        static constexpr std::size_t exct  = 8;
        static constexpr std::size_t exss  = 9;
        static constexpr std::size_t exes  = 10;
        static constexpr std::size_t scr   = 11;
        static constexpr std::size_t gene  = 12;
        static constexpr std::size_t cdssx = 13;
        static constexpr std::size_t cdsex = 14;
        static constexpr std::size_t exfrm = 15;

        static std::initializer_list<std::string> labels()
        {
            return {"bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"};
        }
    };

    class tx_command : public varoom::command
    {
    public:
        tx_command(const string& p_input_filename, const string& p_output_filename)
            : m_input_filename(p_input_filename),
              m_output_filename(p_output_filename)
        {
        }

        using tx_range = std::pair<uint32_t,uint32_t>;
        using tx_tuples = std::vector<refgene::tuple>;

        virtual void operator()()
        {
            std::unordered_map<atom,std::map<tx_range,tx_tuples>> data;

            {
                std::function<void(const refgene::tuple&)> f = [&] (const refgene::tuple& p_entry) mutable {
                    atom ch = std::get<refgene::chr>(p_entry);
                    tx_range x = tx_range(std::get<refgene::txs>(p_entry), std::get<refgene::txe>(p_entry));
                    data[ch][x].push_back(p_entry);
                };
                input_file_holder_ptr globp = files::in(m_input_filename);
                refgene::istream_reader g(**globp, true);
                table_utils::for_each(g, f);
            }

            vector<subtext> exs;
            vector<subtext> exe;
            for (auto itr = data.begin(); itr != data.end(); ++itr)
            {
                ranges_builder R;
                for (auto jtr = itr->second.begin(); jtr != itr->second.end(); ++jtr)
                {
                    ranges_builder::range x(jtr->first.first, jtr->first.second);
                    ranges_builder::range_id x_i = R.insert(x);
                    for (size_t i = 0; i < jtr->second.size(); ++i)
                    {
                        const refgene::tuple& t = jtr->second[i];

                        subtext exs_txt = std::get<refgene::exss>(t);
                        subtext exe_txt = std::get<refgene::exes>(t);
                        exs_txt.split(',', exs);
                        exe_txt.split(',', exe);

                        for (size_t j = 0; j < exs.size(); ++j)
                        {
                            if (exs[j].size() == 0)
                            {
                                continue;
                            }
                            ranges_builder::range ex;
                            ex.first = lexical_cast<size_t>(make_iterator_range(exs[j].first, exs[j].second));
                            ex.second = lexical_cast<size_t>(make_iterator_range(exe[j].first, exe[j].second));
                            ranges_builder::range_id ex_i = R.insert(x);
                            std::cout << itr->first
                                      << '\t' << "<" << ex.first << "," << ex.second << "> -> " << ex_i
                                      << std::endl;
                        }
                    }
                }
            }
        }

    private:
        const string m_input_filename;
        const string m_output_filename;
    };

    class tx_factory : public command_factory
    {
    public:
        tx_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("compile a compact transcript reference");
            opts.add_options()
                ("help,h", "produce help message")
                ("input,i", po::value<string>(), "input filename, defaults to '-' for stdin")
                ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
                ;

            po::positional_options_description pos;
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

            string in_fn = "-";
            if (vm.count("input"))
            {
                in_fn = vm["input"].as<string>();
            }
            params["input"] = in_fn;
            params["output"] = vm["output"].as<string>();

            return params;
        }

        virtual command_ptr create(const json& p_params) const
        {
            if (!p_params.is_object())
            {
                return command_ptr();
            }
            string input_fn = p_params["input"];
            string output_fn = p_params["output"];
            return command_ptr(new tx_command(input_fn, output_fn));
        }
    };
    
    bool reg = command_factory::add("tx", command_factory_ptr(new tx_factory));
}
// namespace anonymous

