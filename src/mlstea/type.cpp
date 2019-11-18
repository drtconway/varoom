#include "varoom/command.hpp"

#include <iostream>
#include <unordered_map>
#include <vector>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>
#include <nlohmann/json.hpp>

#include "varoom/kmers.hpp"
#include "varoom/seq/fasta.hpp"
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

    typedef vector<size_t>  allele_numbers;
    typedef unordered_map<size_t,allele_numbers> locus_alleles;
    typedef unordered_map<kmer,locus_alleles> index_type;

    typedef unordered_map<size_t,size_t> allele_count;
    typedef unordered_map<size_t,allele_count> locus_allele_count;

    class locus_index
    {
    public:
        locus_index(const string& p_locus_file, size_t p_k)
            : m_k(p_k)
        {
            read_locus_file(p_locus_file);
            for (size_t i = 0; i < m_locus_files.size(); ++i)
            {
                vector<string> labels;
                input_file_holder_ptr inp = files::in(m_locus_files[i]);
                for (fasta_reader rx(**inp); rx.more(); ++rx)
                {
                    const fasta_read& r = *rx;
                    size_t n = labels.size();
                    labels.push_back(r.first);
                    vector<kmer> xs;
                    kmers::make(r.second, m_k, xs);
                    for (size_t j = 0; j < xs.size(); ++j)
                    {
                        kmer x = xs[j];
                        allele_numbers& ax = m_index[x][i];
                        if (ax.size() == 0 || ax.back() != n)
                        {
                            ax.push_back(n);
                        }
                    }
                }
                m_locus_labels.push_back(labels);
            }
        }

        size_t lookup(const vector<kmer>& p_xs, locus_allele_count& p_counts) const
        {
            size_t n = 0;
            for (size_t i = 0; i < p_xs.size(); ++i)
            {
                kmer x = p_xs[i];
                auto ix = m_index.find(x);
                if (ix == m_index.end())
                {
                    continue;
                }
                ++n;
                for (auto j = ix->second.begin(); j != ix->second.end(); ++j)
                {
                    allele_count& lx = p_counts[j->first];
                    for (auto k = j->second.begin(); k != j->second.end(); ++k)
                    {
                        lx[*k] += 1;
                    }
                }
            }
            return n;
        }

        json json_counts(const locus_allele_count& p_cx) const
        {
            json r = json::object();
            for (auto i = p_cx.begin(); i != p_cx.end(); ++i)
            {
                size_t l = i->first;
                json s = json::object();
                for (auto j = i->second.begin(); j != i->second.end(); ++j)
                {
                    size_t a = j->first;
                    size_t c = j->second;
                    s[m_locus_labels[l][a]] = c;
                }
                r[m_locus_names[l]] = s;
            }
            return r;
        }

    private:
        void read_locus_file(const string& p_locus_file)
        {
            path loc_p = absolute(p_locus_file);
            path par_p = loc_p.parent_path();

            input_file_holder_ptr inp = files::in(loc_p.string());
            string ln;
            int i = 0;
            while (std::getline(**inp, ln))
            {
                boost::algorithm::trim(ln);
                path l_p(ln);
                path r_p = l_p;
                if (r_p.is_relative())
                {
                    r_p = par_p;
                    r_p /= l_p;
                }
                m_locus_names.push_back(l_p.string());
                m_locus_files.push_back(r_p.string());
                ++i;
            }
        }

        void index_json() const
        {
            json J;
            for (auto i = m_index.begin(); i != m_index.end(); ++i)
            {
                string xS = kmers::render(m_k, i->first);
                for (auto j = i->second.begin(); j != i->second.end(); ++j)
                {
                    string lS = m_locus_names[j->first];
                    J[xS][lS] = range_vector(j->second);
                }
            }
            cout << J << endl;
        }

        static json range_vector(const vector<size_t>& p_nums)
        {
            json r = json::array();
            size_t p = 0;
            size_t n = 0;
            for (size_t i = 0; i < p_nums.size(); ++i)
            {
                size_t x = p_nums[i];
                if (x != p + n)
                {
                    if (n == 1)
                    {
                        r.push_back(p);
                    }
                    else if (n > 1)
                    {
                        json s;
                        s.push_back(p);
                        s.push_back(p+n-1);
                        r.push_back(s);
                    }
                    p = x;
                    n = 0;
                }
                n += 1;
            }
            if (n == 1)
            {
                r.push_back(p);
            }
            else if (n > 1)
            {
                json s;
                s.push_back(p);
                s.push_back(p+n-1);
                r.push_back(s);
            }
            return r;
        }

        const size_t m_k;
        vector<string> m_locus_names;
        vector<string> m_locus_files;
        vector<vector<string>> m_locus_labels;
        index_type m_index;
    };

    class type_command : public varoom::command
    {
    public:
        type_command(const strings& p_input_filenames,
                            const string& p_output_filename,
                            const string& p_locus_file)
            : m_input_filenames(p_input_filenames),
              m_output_filename(p_output_filename),
              m_locus_file(p_locus_file)
        {
        }

        virtual void operator()()
        {
            const size_t K = 21;
            locus_index L(m_locus_file, K);

            locus_allele_count cx;

            vector<kmer> lhsFwd;
            vector<kmer> lhsRev;
            vector<kmer> rhsFwd;
            vector<kmer> rhsRev;
            for (size_t i = 0; i + 1 < m_input_filenames.size(); ++i)
            {
                input_file_holder_ptr inp1 = files::in(m_input_filenames[i]);
                input_file_holder_ptr inp2 = files::in(m_input_filenames[i+1]);
                fastq_reader r1(**inp1);
                fastq_reader r2(**inp2);
                for (fastq_pair_reader r1r2(r1, r2); r1r2.more(); ++r1r2)
                {
                    kmers::make(std::get<1>((*r1r2).first), K, lhsFwd, lhsRev);
                    kmers::make(std::get<1>((*r1r2).second), K, rhsFwd, rhsRev);

                    L.lookup(lhsFwd, cx);
                    L.lookup(rhsRev, cx);

                    L.lookup(rhsFwd, cx);
                    L.lookup(lhsRev, cx);

                }
            }
            cout << L.json_counts(cx) << endl;
        }

    private:

        const strings m_input_filenames;
        const string m_output_filename;
        const string m_locus_file;
    };

    class type_factory : public command_factory
    {
    public:
        type_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const vector<string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("count bases");
            opts.add_options()
                ("help,h", "produce help message")
                ("loci,l", po::value<string>(), "only include counts from regions in the named BED file")
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

            if (vm.count("loci"))
            {
                params["loci"] = vm["loci"].as<string>();
            }

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
            string loci_fn = p_params["loci"];
            return command_ptr(new type_command(inputs_fn, output_fn, loci_fn));
        }
    };
    
    bool reg = command_factory::add("type", command_factory_ptr(new type_factory));
}
// namespace anonymous

