#include "varoom/command.hpp"

#include <boost/format.hpp>

#include "varoom/kmers.hpp"
#include "varoom/seq/fasta.hpp"
#include "varoom/seq/fastq.hpp"
#include "varoom/seq/fastq_pair.hpp"
#include "varoom/util/codec8.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/murmur3.hpp"
#include "varoom/util/subtext.hpp"

#include <nlohmann/json.hpp>

#include <queue>
#include <unordered_map>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::seq;

namespace // anonymous
{
    using json = nlohmann::json;
    typedef fastq_pair<fastq_reader,fastq_reader> fastq_pair_reader;
    typedef fastq_pair_reader::item_type read_pair;

    class byte_iterator
    {
    public:
        byte_iterator(std::istream& p_in)
            : m_in(p_in)
        {
            m_valid = next(m_x);
        }

        bool valid() const
        {
            return m_valid;
        }

        uint8_t operator*() const
        {
            return m_x;
        }

        byte_iterator& operator++()
        {
            m_valid = next(m_x);
            return *this;
        }

        bool next(uint8_t& p_x)
        {
            char c;
            if (m_in.get(c))
            {
                p_x = static_cast<uint8_t>(c);
                return true;
            }
            return false;
        }

    private:
        std::istream& m_in;
        bool m_valid;
        uint8_t m_x;
    };

    class kmer_istream
    {
    public:
        kmer_istream(const size_t& p_label, const std::string& p_filename)
            : m_label(p_label), m_file_holder(files::in(p_filename)), m_bytes(**m_file_holder), m_x(0)
        {
            m_valid = m_bytes.valid();
            if (m_valid)
            {
                m_x += util::codec8::decode(m_bytes);
            }
        }

        bool valid() const
        {
            return m_valid;
        }

        const size_t& label() const
        {
            return m_label;
        }

        kmer operator*() const
        {
            return m_x;
        }

        kmer_istream& operator++()
        {
            if (m_bytes.valid())
            {
                m_x += util::codec8::decode(m_bytes);
            }
            else
            {
                m_valid = false;
            }
            return *this;
        }

        bool operator<(const kmer_istream& other) const
        {
            if (!valid())
            {
                return false;
            }
            if (!other.valid())
            {
                return true;
            }
            return *(*this) > *other;
        }

    private:
        const size_t m_label;
        input_file_holder_ptr m_file_holder;
        byte_iterator m_bytes;
        kmer m_x;
        bool m_valid;
    };
    using kmer_stream_ptr = std::shared_ptr<kmer_istream>;

    class vbyte_istream
    {
    public:
        vbyte_istream(const std::string& p_filename)
            : m_file_holder(files::in(p_filename)), m_bytes(**m_file_holder), m_x(0)
        {
            m_valid = m_bytes.valid();
            if (m_valid)
            {
                m_x += util::codec8::decode(m_bytes);
            }
        }

        bool valid() const
        {
            return m_valid;
        }

        uint64_t operator*() const
        {
            return m_x;
        }

        vbyte_istream& operator++()
        {
            if (m_bytes.valid())
            {
                m_x = util::codec8::decode(m_bytes);
            }
            else
            {
                m_valid = false;
            }
            return *this;
        }

        bool operator<(const vbyte_istream& other) const
        {
            if (!valid())
            {
                return false;
            }
            if (!other.valid())
            {
                return true;
            }
            return *(*this) > *other;
        }

    private:
        input_file_holder_ptr m_file_holder;
        byte_iterator m_bytes;
        kmer m_x;
        bool m_valid;
    };
    using kmer_stream_ptr = std::shared_ptr<kmer_istream>;

    class kmer_ostream
    {
    public:
        kmer_ostream(const std::string& p_filename)
            : m_file_holder(files::out(p_filename)), m_out(**m_file_holder), m_x(0)
        {
        }

        kmer_ostream& push_back(const kmer& p_x)
        {
            kmer d = p_x - m_x;
            m_bs.clear();
            util::codec8::encode(d, m_bs);
            m_out.write(reinterpret_cast<const char*>(&m_bs[0]), m_bs.size());
            m_x = p_x;
            return *this;
        }

    private:
        output_file_holder_ptr m_file_holder;
        std::ostream& m_out;
        kmer m_x;
        std::vector<uint8_t> m_bs;
    };

    class vbyte_ostream
    {
    public:
        vbyte_ostream(const std::string& p_filename)
            : m_file_holder(files::out(p_filename)), m_out(**m_file_holder)
        {
        }

        vbyte_ostream& push_back(const kmer& p_x)
        {
            m_bs.clear();
            util::codec8::encode(p_x, m_bs);
            m_out.write(reinterpret_cast<const char*>(&m_bs[0]), m_bs.size());
            return *this;
        }

    private:
        output_file_holder_ptr m_file_holder;
        std::ostream& m_out;
        std::vector<uint8_t> m_bs;
    };

    class kmerize_command : public varoom::command
    {
    public:
        kmerize_command(const size_t& p_k, const string& p_input_filename, const string& p_output_filename)
            : m_k(p_k), m_input_filename(p_input_filename), m_output_filename(p_output_filename)
        {
        }

        void create_kmers()
        {
            vector<kmer> xs;
            input_file_holder_ptr inp = files::in(m_input_filename);
            size_t n = 0;
            map<string,string> names;
            for (fasta_reader rx(**inp); rx.more(); ++rx, ++n)
            {
                const fasta_read& r = *rx;
                subtext tx(r.first);
                vector<subtext> txs;
                tx.split(' ', txs);
                string nm = string(txs[0]);
                std::cerr << "got sequence: " << nm << std::endl;
                xs.reserve(r.second.size());
                kmers::make(r.second, m_k, xs);

                std::sort(xs.begin(), xs.end());
                xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
                
                std::string fn = str(format("%s-%03d.vbyte") % m_output_filename % n);
                names[nm] = fn;
                output_file_holder_ptr outp = files::out(fn);
                std::ostream& out = **outp;

                kmer px = 0;
                std::vector<uint8_t> bs;
                for (size_t i = 0; i < xs.size(); ++i)
                {
                    kmer x = xs[i];
                    kmer d = x - px;
                    bs.clear();
                    util::codec8::encode(d, bs);
                    out.write(reinterpret_cast<const char*>(&bs[0]), bs.size());
                    px = x;
                }
            }
            {
                std::string fn = str(format("%s-toc.json") % m_output_filename);
                output_file_holder_ptr outp = files::out(fn);
                std::ostream& out = **outp;
                out << json(names);
            }
        }

        void merge_kmers()
        {
            map<string,string> names;
            {
                std::string fn = str(format("%s-toc.json") % m_output_filename);
                input_file_holder_ptr inp = files::in(fn);
                std::istream& in = **inp;
                json j;
                in >> j;
                for (auto itr = j.begin(); itr != j.end(); ++itr)
                {
                    names[itr.key()] = itr.value();
                }
            }

            auto cmp = [](const kmer_stream_ptr& p_lhs, const kmer_stream_ptr& p_rhs) {
                return *p_lhs < *p_rhs;
            };
            priority_queue<kmer_stream_ptr, vector<kmer_stream_ptr>, decltype(cmp)> Q(cmp);
            size_t n = 0;
            vector<string> labels;
            for (auto itr = names.begin(); itr != names.end(); ++itr, ++n)
            {
                labels.push_back(itr->first);
                kmer_stream_ptr ksp(new kmer_istream(n,  itr->second));
                if (ksp->valid())
                {
                    Q.push(ksp);
                }
            }

            kmer_ostream kout(m_output_filename + "-merged.vbyte");
            vbyte_ostream vout(m_output_filename + "-chroms.vbyte");

            vector<size_t> cur;
            while (Q.size() > 0)
            {
                cur.clear();
                kmer_stream_ptr ksp = Q.top();
                Q.pop();
                kmer x = **ksp;
                cur.push_back(ksp->label());
                ++(*ksp);
                if (ksp->valid())
                {
                    Q.push(ksp);
                }
                while (Q.size() > 0 && **Q.top() == x)
                {
                    kmer_stream_ptr jsp = Q.top();
                    Q.pop();
                    cur.push_back(jsp->label());
                    ++(*jsp);
                    if (jsp->valid())
                    {
                        Q.push(jsp);
                    }
                }
                size_t w = (cur.size() == 1 ? cur[0] + 1 : 0);
                kout.push_back(x);
                vout.push_back(w);
            }
        }

        void classify_reads()
        {
            unordered_map<kmer,vector<size_t>> idx1;
            unordered_map<kmer,vector<size_t>> idx2;

            vector<kmer> lhsFwd;
            vector<kmer> lhsRev;
            vector<kmer> rhsFwd;
            vector<kmer> rhsRev;
            if (1)
            {
                input_file_holder_ptr inp1 = files::in("/data/work/her2-sra/ERR863738_1.fastq.gz");
                input_file_holder_ptr inp2 = files::in("/data/work/her2-sra/ERR863738_1.fastq.gz");
                fastq_reader r1(**inp1);
                fastq_reader r2(**inp2);
                size_t n = 0;
                for (fastq_pair_reader r1r2(r1, r2); r1r2.more(); ++r1r2, ++n)
                {
                    kmers::make(std::get<1>((*r1r2).first), m_k, lhsFwd, lhsRev);
                    kmers::make(std::get<1>((*r1r2).second), m_k, rhsFwd, rhsRev);

                    lhsFwd.insert(lhsFwd.end(), rhsRev.begin(), rhsRev.end());
                    std::sort(lhsFwd.begin(), lhsFwd.end());
                    lhsFwd.erase(std::unique(lhsFwd.begin(), lhsFwd.end()), lhsFwd.end());
                    for (auto itr = lhsFwd.begin(); itr != lhsFwd.end(); ++itr)
                    {
                        idx1[*itr].push_back(n);
                    }

                    rhsFwd.insert(rhsFwd.end(), lhsRev.begin(), lhsRev.end());
                    std::sort(rhsFwd.begin(), rhsFwd.end());
                    rhsFwd.erase(std::unique(rhsFwd.begin(), rhsFwd.end()), rhsFwd.end());
                    for (auto itr = rhsFwd.begin(); itr != rhsFwd.end(); ++itr)
                    {
                        idx2[*itr].push_back(n);
                    }

                    constexpr size_t N = (1ULL << 16) - 1;
                    if ((n & N) == N)
                    {
                        std::cerr << "scanning!" << std::endl;

                        kmer_istream kin(0, m_output_filename + "-merged.vbyte");
                        vbyte_istream vin(m_output_filename + "-chroms.vbyte");

                        vector<set<uint64_t>> xs1(N+1);
                        vector<set<uint64_t>> xs2(N+1);
            
                        constexpr size_t J = 996011*13;
                        size_t j = J;
                        while (kin.valid())
                        {
                            kmer x = *kin;
                            ++kin;
                            {
                                --j;
                                if (j == 0)
                                {
                                    std::cerr << kmers::render(m_k, x) << std::endl;
                                    j += J;
                                }
                            }
                            uint64_t c = *vin;
                            ++vin;

                            auto f1 = idx1.find(x);
                            if (f1 != idx1.end())
                            {
                                for (auto itr = f1->second.begin(); itr != f1->second.end(); ++itr)
                                {
                                    xs1[*itr & N].insert(c);
                                }
                            }

                            auto f2 = idx2.find(x);
                            if (f2 != idx2.end())
                            {
                                for (auto itr = f2->second.begin(); itr != f2->second.end(); ++itr)
                                {
                                    xs2[*itr & N].insert(c);
                                }
                            }
                        }

                        for (size_t i = 0; i < xs1.size(); ++i)
                        {
                            if (xs1[i].size() > 0 || xs2[i].size() > 0)
                            {
                                vector<uint64_t> ys1(xs1[i].begin(), xs1[i].end());
                                vector<uint64_t> ys2(xs2[i].begin(), xs2[i].end());
                                std::cout << i << '\t' << json(ys1) << '\t' << json(ys2) << endl;
                            }
                        }

                        idx1.clear();
                        idx2.clear();
                        break;
                    }
                }
            }
        }

        virtual void operator()()
        {
            classify_reads();
            //merge_kmers();
        }

    private:
        const size_t m_k;
        const string m_input_filename;
        const string m_output_filename;
    };

    class kmerize_factory : public command_factory
    {
    public:
        kmerize_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("compile a compact genome reference");
            opts.add_options()
                ("help,h", "produce help message")
                ("input,i", po::value<string>()->default_value("-"), "input filename, defaults to '-' for stdin")
                ("kmer,k", po::value<size_t>()->default_value(25), "kmer size")
                ("output,o", po::value<string>(), "output filename")
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

            params["kmer"] = vm["kmer"].as<size_t>();
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
            size_t k = p_params["kmer"];
            string input_fn = p_params["input"];
            string output_fn = p_params["output"];
            return command_ptr(new kmerize_command(k, input_fn, output_fn));
        }
    };
    
    bool reg = command_factory::add("kmerize", command_factory_ptr(new kmerize_factory));
}
// namespace anonymous


