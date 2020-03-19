#include "varoom/command.hpp"

#include <map>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <boost/format.hpp>

#include "varoom/kmers.hpp"
#include "varoom/seq/fasta.hpp"
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
    string human(uint64_t x)
    {
        if (x < 1000)
        {
            return str(format("%d") % x);
        }
        double y = x / 1000.0;
        if (y < 1000)
        {
            return str(format("%2.2fK") % y);
        }
        y = y / 1000.0;
        if (y < 1000)
        {
            return str(format("%2.2fM") % y);
        }
        y = y / 1000.0;
        if (y < 1000)
        {
            return str(format("%2.2fG") % y);
        }
        y = y / 1000.0;
        return str(format("%2.2fT") % y);
    }

    using strings = vector<string>;
    using fastq_pair_reader = fastq_pair<fastq_reader,fastq_reader>;
    using read_pair = fastq_pair_reader::item_type;

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

    bool hits(const vector<kmer>& p_idx, vector<kmer>& p_xs)
    {
        std::sort(p_xs.begin(), p_xs.end());
        auto lb = p_idx.begin();
        auto ub = p_idx.end();
        for (size_t i = 0; i < p_xs.size(); ++i)
        {
            kmer x = p_xs[i];
            auto itr = std::lower_bound(lb, ub, x);
            if (itr == ub)
            {
                return false;
            }
            if (*itr == x)
            {
                return true;
            }
            lb = itr;
        }
        return false;
    }

    class pulldown_command : public varoom::command
    {
    public:
        pulldown_command(const string& p_reference_filename, const strings& p_input_filenames, const strings& p_output_filenames)
            : m_reference_filename(p_reference_filename),
              m_input_filenames(p_input_filenames),
              m_output_filenames(p_output_filenames)
        {
        }

        virtual void operator()()
        {
            constexpr size_t K = 25;
            constexpr size_t J = 4;
            constexpr size_t Q = 1;

            vector<kmer> R;
            {
                std::cerr << "reading " << m_reference_filename << std::endl;
                unordered_set<kmer> R0;
                input_file_holder_ptr inp = files::in(m_reference_filename);
                vector<kmer> xs;
                for (fasta_reader rx(**inp); rx.more(); ++rx)
                {
                    const fasta_read& r = *rx;
                    xs.clear();
                    kmers::make(r.second, K, xs);
                    for (size_t i = 0; i < xs.size(); ++i)
                    {
                        kmer x = xs[i];
                        if (1)
                        {
                            uint32_t a = murmur3(17).update(x)();
                            if ((a % J) >= Q)
                            {
                                // skip this kmer.
                                continue;
                            }
                        }
                        R0.insert(x);
                    }
                }
                std::cerr << "done." << std::endl;
                R.insert(R.end(), R0.begin(), R0.end());
            }
            std::cerr << "sorting..." << std::endl;
            std::sort(R.begin(), R.end());
            std::cerr << "done." << std::endl;


            size_t rn = 0;
            size_t lhsFwdCnt = 0;
            size_t lhsRevCnt = 0;
            size_t rhsFwdCnt = 0;
            size_t rhsRevCnt = 0;

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
                for (fastq_pair_reader r1r2(r1, r2); r1r2.more(); ++r1r2, ++rn)
                {
                    if (1)
                    {
                        uint32_t a = murmur3(17).update(rn)();
                        if ((a % 1024) < 1)
                        {
                            std::cerr << human(rn) << '\r';
                            std::cerr.flush();
                        }
                    }
                    
                    kmers::make(std::get<1>((*r1r2).first), K, lhsFwd, lhsRev);
                    kmers::make(std::get<1>((*r1r2).second), K, rhsFwd, rhsRev);

                    if (hits(R, lhsFwd))
                    {
                        lhsFwdCnt += 1;
                    }
                    if (hits(R, lhsRev))
                    {
                        lhsRevCnt += 1;
                    }
                    if (hits(R, rhsFwd))
                    {
                        rhsFwdCnt += 1;
                    }
                    if (hits(R, rhsRev))
                    {
                        rhsRevCnt += 1;
                    }
                }
            }

            std::cout << rn
                      << '\t' << lhsFwdCnt
                      << '\t' << lhsRevCnt
                      << '\t' << rhsFwdCnt
                      << '\t' << rhsRevCnt
                      << std::endl;
        }

    private:
        const string m_reference_filename;
        const strings m_input_filenames;
        const strings m_output_filenames;
    };

    class pulldown_factory : public command_factory
    {
    public:
        pulldown_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("compile a compact genome reference");
            opts.add_options()
                ("help,h", "produce help message")
                ("reference,r", po::value<string>(), "reference sequence filename")
                ("input,i", po::value<strings>(), "input filenames")
                ("output,o", po::value<strings>(), "output filename, defaults to '-' for stdout")
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

            params["reference"] = vm["reference"].as<string>();
            params["input"] = vm["input"].as<strings>();
            params["output"] = vm["output"].as<strings>();

            return params;
        }

        virtual command_ptr create(const json& p_params) const
        {
            if (!p_params.is_object())
            {
                return command_ptr();
            }
            string ref_fn = p_params["reference"];
            strings inputs_fn = p_params["input"];
            strings outputs_fn = p_params["output"];
            return command_ptr(new pulldown_command(ref_fn, inputs_fn, outputs_fn));
        }
    };
    
    bool reg = command_factory::add("pulldown", command_factory_ptr(new pulldown_factory));
}
// namespace anonymous

