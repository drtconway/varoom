#include "varoom/command.hpp"

#include <map>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <iostream>

#include "varoom/kmers.hpp"
#include "varoom/seq/fastq.hpp"
#include "varoom/seq/fastq_pair.hpp"
#include "varoom/util/murmur3.hpp"
#include "varoom/util/files.hpp"

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::seq;

namespace // anonymous
{
    class minhash
    {
    public:
        struct hash_comparitor
        {
            bool operator()(const kmer& p_x, const kmer& p_y) const
            {
                uint32_t a = murmur3(19).update(p_x)();
                uint32_t b = murmur3(19).update(p_y)();
                return a < b;
            }
        };

        minhash(size_t p_max_size)
            : m_max_size(p_max_size)
        {
        }

        minhash& add(const kmer& p_x)
        {
            if (m_items.count(p_x) > 0)
            {
                return *this;
            }

            m_items.insert(p_x);
            m_queue.push(p_x);

            while (m_items.size() > m_max_size)
            {
                kmer x = m_queue.top();
                m_queue.pop();
                m_items.erase(x);
            }

            return *this;
        }

        const std::unordered_set<kmer>& items() const
        {
            return m_items;
        }

    private:
        const size_t m_max_size;
        std::unordered_set<kmer> m_items;
        std::priority_queue<kmer, std::vector<kmer>, hash_comparitor> m_queue;
    };

    class sketch_command : public varoom::command
    {
    public:
        sketch_command(const uint64_t& p_size, const string& p_input_filename, const string& p_output_filename)
            : m_size(p_size),
              m_input_filename(p_input_filename),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            input_file_holder_ptr inp1 = files::in(m_input_filename);

            constexpr size_t K = 25;
            constexpr size_t J = 16;
            constexpr size_t Q = 1;

            vector<kmer> fwd;
            vector<kmer> rev;
            unordered_map<kmer,size_t> X;
            for (fastq_reader r1(**inp1); r1.more(); ++r1)
            {
                kmers::make(std::get<1>(*r1), K, fwd, rev);
                for (size_t i = 0; i < fwd.size(); ++i)
                {
                    kmer x = std::min(fwd[i], rev[i]);
                    uint32_t a = murmur3(17).update(x)();
                    if ((a % J) >= Q)
                    {
                        // skip this kmer.
                        continue;
                    }
                    X[x] += 1;
                }
            }
            map<size_t,size_t> H;
            for (auto itr = X.begin(); itr != X.end(); ++itr)
            {
                H[itr->second] += 1;
            }
            for (auto itr = H.begin(); itr != H.end(); ++itr)
            {
                std::cout << itr->first << '\t' << itr->second << std::endl;
            }
        }

    private:
        const uint64_t m_size;
        const string m_input_filename;
        const string m_output_filename;
    };

    class sketch_factory : public command_factory
    {
    public:
        sketch_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("compile a compact genome reference");
            opts.add_options()
                ("help,h", "produce help message")
                ("size,z", po::value<uint64_t>()->default_value(65536), "maximum number of k-mers to keep [default: 65536]")
                ("input,i", po::value<string>()->default_value("-"), "input filename, defaults to '-' for stdin")
                ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
                ;

            po::positional_options_description pos;

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

            params["size"] = vm["size"].as<uint64_t>();
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
            uint64_t z = p_params["size"];
            string input_fn = p_params["input"];
            string output_fn = p_params["output"];
            return command_ptr(new sketch_command(z, input_fn, output_fn));
        }
    };
    
    bool reg = command_factory::add("sketch", command_factory_ptr(new sketch_factory));
}
// namespace anonymous

