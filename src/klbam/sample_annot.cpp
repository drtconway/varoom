#include "varoom/command.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/text.hpp"
#include "varoom/util/table.hpp"
#include "varoom/seq/index.hpp"
#include "varoom/vcf/vcf_parser.hpp"
#include "varoom/klbam/types.hpp"

#include <cmath>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::vcf;
using namespace varoom::klbam;

namespace // anonymous
{
    struct annot : table<chrom,uint32_t,
                         uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,uint32_t,string,
                         double,double,
                         uint32_t,double,double,
                         string,string,string>
    {
        using tuple = tuple_type;
        using table = basic_inmemory_table;
        using istream_reader = basic_istream_table;
        using ostream_writer = basic_ostream_table;
        using read_iterator = basic_read_iterator;
        using write_iterator = basic_write_iterator;

        static constexpr std::size_t chr = 0;
        static constexpr std::size_t pos = 1;
        static constexpr std::size_t nA  = 2;
        static constexpr std::size_t nC  = 3;
        static constexpr std::size_t nG  = 4;
        static constexpr std::size_t nT  = 5;
        static constexpr std::size_t nN  = 6;
        static constexpr std::size_t nX0 = 7;
        static constexpr std::size_t nX1 = 8;
        static constexpr std::size_t nI  = 9;
        static constexpr std::size_t nD  = 10;
        static constexpr std::size_t oth = 11;
        static constexpr std::size_t kld  = 12;
        static constexpr std::size_t pval = 13;
        static constexpr std::size_t cov = 14;
        static constexpr std::size_t raf = 15;
        static constexpr std::size_t vaf = 16;
        static constexpr std::size_t ref = 17;
        static constexpr std::size_t alt = 18;
        static constexpr std::size_t exac = 19;

        static std::initializer_list<std::string> labels()
        {
            return {"chr", "pos", "nA", "nC", "nG", "nT", "nN", "nX0", "nX1", "nI", "nD", "indels",
                    "kld", "pval", "cov", "raf", "vaf", "ref", "alt", "exac"};
        }
    };

    class annotator : public vcf::vcf_handler
    {
    public:
        annotator(map<locus,annot::tuple>& p_items)
            : m_items(p_items)
        {
        }

        virtual void operator()(const std::string& p_chr,
                                const std::int64_t& p_pos,
                                const std::string& p_id,
                                const std::string& p_ref,
                                const std::string& p_alt,
                                const double& p_qual,
                                const std::string& p_filter,
                                const lazy<vcf_info>& p_info,
                                const lazy<std::vector<vcf_info>>& p_genotypes)
        {
            locus l(p_chr, p_pos);
            auto i = m_items.find(l);
            if (i == m_items.end())
            {
                return;
            }
            annot::tuple& t = i->second;
            const vcf_info& ifo = p_info.get();
            if (ifo.contains("AF"))
            {
                std::get<annot::exac>(t) = ifo["AF"];
            }
            if (p_ref[0] == 'A')
            {
                std::get<annot::raf>(t) = double(std::get<annot::nA>(t)) / double(std::get<annot::cov>(t));
                std::get<annot::vaf>(t) = double(std::max({
                        std::get<annot::nC>(t), std::get<annot::nG>(t), std::get<annot::nT>(t),
                        std::get<annot::nI>(t), std::get<annot::nD>(t)}))
                     / double(std::get<annot::cov>(t));
            }
            else if (p_ref[0] == 'C')
            {
                std::get<annot::raf>(t) = double(std::get<annot::nC>(t)) / double(std::get<annot::cov>(t));
                std::get<annot::vaf>(t) = double(std::max({
                        std::get<annot::nA>(t), std::get<annot::nG>(t), std::get<annot::nT>(t),
                        std::get<annot::nI>(t), std::get<annot::nD>(t)}))
                     / double(std::get<annot::cov>(t));
            }
            else if (p_ref[0] == 'G')
            {
                std::get<annot::raf>(t) = double(std::get<annot::nG>(t)) / double(std::get<annot::cov>(t));
                std::get<annot::vaf>(t) = double(std::max({
                        std::get<annot::nA>(t), std::get<annot::nC>(t), std::get<annot::nT>(t),
                        std::get<annot::nI>(t), std::get<annot::nD>(t)}))
                     / double(std::get<annot::cov>(t));
            }
            else if (p_ref[0] == 'T')
            {
                std::get<annot::raf>(t) = double(std::get<annot::nT>(t)) / double(std::get<annot::cov>(t));
                std::get<annot::vaf>(t) = double(std::max({
                        std::get<annot::nA>(t), std::get<annot::nC>(t), std::get<annot::nG>(t),
                        std::get<annot::nI>(t), std::get<annot::nD>(t)}))
                     / double(std::get<annot::cov>(t));
            }
        }

    private:
        map<locus,annot::tuple>& m_items;
    };

    class sample_annot_command : public varoom::command
    {
    public:
        sample_annot_command(const string& p_counts_filename,
                            const string& p_scores_filename,
                            const string& p_genome_filename,
                            const string& p_output_filename)
            : m_counts_filename(p_counts_filename),
              m_scores_filename(p_scores_filename),
              m_genome_filename(p_genome_filename),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            varoom::seq::index I(m_genome_filename, 20);

            map<locus,annot::tuple> items;
            {
                std::function<void(const counts::tuple&,const scores::tuple&)> f =
                    [&] (const counts::tuple& p_counts,const scores::tuple& p_score) mutable
                {
                    locus l(std::get<counts::chr>(p_counts), std::get<counts::pos>(p_counts));
                    if (!(std::get<scores::chr>(p_score) == l.first && std::get<scores::pos>(p_score) == l.second))
                    {
                        throw std::runtime_error("counts and scores out of sync!");
                    }
                    annot::tuple t;
                    table_utils::copy<0,12>(p_counts, t);
                    table_utils::copy<2,4,10>(p_score, t);
                    std::get<annot::cov>(t) = std::get<annot::nA>(t) + std::get<annot::nC>(t)
                                              + std::get<annot::nG>(t) + std::get<annot::nT>(t)
                                              + std::get<annot::nI>(t) + std::get<annot::nD>(t);

                    std::get<annot::ref>(t) = I.get(l.first, l.second-1, l.second);

                    vector<pair<uint32_t,string>> cts;
                    if (std::get<annot::ref>(t) != "A" && std::get<annot::ref>(t) != "a")
                    {
                        cts.push_back(make_pair(std::get<annot::nA>(t), "A"));
                    }
                    else
                    {
                        std::get<annot::raf>(t) = (double)std::get<annot::nA>(t) / (double)std::get<annot::cov>(t);
                    }
                    if (std::get<annot::ref>(t) != "C" && std::get<annot::ref>(t) != "c")
                    {
                        cts.push_back(make_pair(std::get<annot::nC>(t), "C"));
                    }
                    else
                    {
                        std::get<annot::raf>(t) = (double)std::get<annot::nC>(t) / (double)std::get<annot::cov>(t);
                    }
                    if (std::get<annot::ref>(t) != "G" && std::get<annot::ref>(t) != "g")
                    {
                        cts.push_back(make_pair(std::get<annot::nG>(t), "G"));
                    }
                    else
                    {
                        std::get<annot::raf>(t) = (double)std::get<annot::nG>(t) / (double)std::get<annot::cov>(t);
                    }
                    if (std::get<annot::ref>(t) != "T" && std::get<annot::ref>(t) != "t")
                    {
                        cts.push_back(make_pair(std::get<annot::nT>(t), "T"));
                    }
                    else
                    {
                        std::get<annot::raf>(t) = (double)std::get<annot::nT>(t) / (double)std::get<annot::cov>(t);
                    }
                    cts.push_back(make_pair(std::get<annot::nI>(t), "I"));
                    cts.push_back(make_pair(std::get<annot::nD>(t), "D"));
                    std::sort(cts.rbegin(), cts.rend());

                    std::get<annot::alt>(t) = cts.begin()->second;
                    std::get<annot::vaf>(t) = (double)(cts.begin()->first) / (double)std::get<annot::cov>(t);
                    
                    //std::cerr << l.first << '\t' << l.second << '\t' << std::get<annot::ref>(t) << std::endl;

                    items[l] = t;
                };
                input_file_holder_ptr cinp = files::in(m_counts_filename);
                counts::istream_reader cs(**cinp, true);
                input_file_holder_ptr sinp = files::in(m_scores_filename);
                scores::istream_reader ss(**sinp, true);
                table_utils::for_each_2(cs, ss, f);
            }

            output_file_holder_ptr outp = files::out(m_output_filename);
            annot::ostream_writer out(**outp, annot::labels());
            for (auto itr = items.begin(); itr != items.end(); ++itr)
            {
                const annot::tuple& t = itr->second;
                if (std::get<annot::pval>(t) > 1e-3)
                {
                    continue;
                }
                //std::cerr << std::get<annot::pval>(t) << '\t' << std::get<annot::ref>(t).size() << std::endl;
                out << t;
            }
        }

    private:
        const string m_counts_filename;
        const string m_scores_filename;
        const string m_genome_filename;
        const string m_output_filename;

    };

    class sample_annot_factory : public command_factory
    {
    public:
        sample_annot_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("annotate results");
            opts.add_options()
                ("help,h", "produce help message")
                ("genome,g", po::value<string>(), "filename for genome table-of-contents")
                ("counts,c", po::value<string>(), "input counts file")
                ("scores,s", po::value<string>(), "input scores file")
                ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
                ;

            po::positional_options_description pos;
            pos.add("genome", 1);
            pos.add("counts", 1);
            pos.add("scores", 1);
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

            params["counts"] = vm["counts"].as<string>();
            params["scores"] = vm["scores"].as<string>();
            params["genome"] = vm["genome"].as<string>();
            params["output"] = vm["output"].as<string>();

            return params;
        }

        virtual command_ptr create(const json& p_params) const
        {
            if (!p_params.is_object())
            {
                return command_ptr();
            }
            string counts_fn = p_params["counts"];
            string scores_fn = p_params["scores"];
            string genome = p_params["genome"];
            string output_fn = p_params["output"];
            return command_ptr(new sample_annot_command(counts_fn, scores_fn, genome, output_fn));
        }
    };
    
    bool reg = command_factory::add("annot", command_factory_ptr(new sample_annot_factory));
}

// namespace anonymous
