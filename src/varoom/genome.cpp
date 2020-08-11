#include "varoom/command.hpp"
#include "varoom/seq/fasta.hpp"
#include "varoom/util/files.hpp"
#include "varoom/kmers.hpp"

#include <unordered_map>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::seq;

namespace // anonymous
{

    using strings = vector<string>;
    using count_and_lastpos = std::pair<size_t,size_t>;

    template <typename K, typename V, typename U>
    class buffered_unordered_map : public unordered_map<K,V>
    {
    public:
        using item_type = std::pair<K, U>;
        using items_type = std::vector<item_type>;
        using map_type = std::unordered_map<K,V>;
        using acceptor_type = std::function<void(map_type&, typename items_type::const_iterator, typename items_type::const_iterator)>;
        static constexpr size_t BufferMin = (1ULL << 16);

        buffered_unordered_map(acceptor_type p_acceptor)
            : m_buffer_max(BufferMin), m_acceptor(p_acceptor)
        {
        }

        buffered_unordered_map& push_back(const item_type& p_item)
        {
            m_buffer.push_back(p_item);
            if (m_buffer.size() >= m_buffer_max)
            {
                flush();
            }
            return *this;
        }

        void flush()
        {
            std::sort(m_buffer.begin(), m_buffer.end());
            auto fst = m_buffer.begin();
            auto lst = m_buffer.begin();
            for (auto jtr = m_buffer.begin(); jtr != m_buffer.end(); ++jtr)
            {
                if (lst->first != jtr->first)
                {
                    m_acceptor(*this, fst, jtr);
                    fst = jtr;
                }
                lst = jtr;
            }
            if (fst != m_buffer.end())
            {
                m_acceptor(*this, fst, m_buffer.end());
            }
            m_buffer.clear();
        }

    private:
        items_type m_buffer;
        size_t m_buffer_max;
        acceptor_type m_acceptor;
    };

    void count_non_overlapping(unordered_map<kmer,size_t>& p_map,
                               vector<pair<kmer,size_t>>::const_iterator p_begin, vector<pair<kmer,size_t>>::const_iterator p_end)
    {
        const size_t K = 20;
        size_t last = 0;
        for (auto itr = p_begin; itr != p_end; ++itr)
        {
            if (itr->second < last + K)
            {
                continue;
            }
            p_map[itr->first] += 1;
            last = itr->second;
        }
    }

    class genome_command : public varoom::command
    {
    public:
        genome_command(const strings& p_input_filenames, const string& p_output_filename)
            : m_input_filenames(p_input_filenames),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            // varoom::seq::index::make(m_input_filenames, m_output_filename);
            const size_t K = 20;
            buffered_unordered_map<kmer,size_t, size_t> symCounts(count_non_overlapping);
            for (size_t i = 0; i < m_input_filenames.size(); ++i)
            {
                std::cerr << "processing " << m_input_filenames[i] << std::endl;
                input_file_holder_ptr inp = files::in(m_input_filenames[i]);
                for (varoom::seq::fasta_reader r(**inp); r.more(); ++r)
                {
                    const fasta_read rd = *r;
                    kmers::make(rd.second, K, [&](kmer_and_pos xp) mutable {
                        symCounts.push_back(xp);
                    });
                }
            }
            symCounts.flush();
            for (auto itr = symCounts.begin(); itr != symCounts.end(); ++itr)
            {
                if (itr->second > 1)
                {
                    cout << kmers::render(K, itr->first) << '\t' << itr->second << endl;
                }
            }
        }

    private:
        const strings m_input_filenames;
        const string m_output_filename;
    };

    class genome_factory : public command_factory
    {
    public:
        genome_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("compile a compact genome reference");
            opts.add_options()
                ("help,h", "produce help message")
                ("input,i", po::value<strings>(), "input filename, defaults to '-' for stdin")
                ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
                ;

            po::positional_options_description pos;
            pos.add("output", 1);
            pos.add("input", -1);

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

            strings ss;
            if (vm.count("input"))
            {
                ss = vm["input"].as<strings>();
            }
            params["input"] = ss;
            params["output"] = vm["output"].as<string>();

            return params;
        }

        virtual command_ptr create(const json& p_params) const
        {
            if (!p_params.is_object())
            {
                return command_ptr();
            }
            strings input_fn = p_params["input"];
            string output_fn = p_params["output"];
            return command_ptr(new genome_command(input_fn, output_fn));
        }
    };
    
    bool reg = command_factory::add("genome", command_factory_ptr(new genome_factory));
}
// namespace anonymous

