#include "varoom/command.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/tsv.hpp"
#include "varoom/util/typed_tsv.hpp"

#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;
using namespace varoom;

namespace // anonymous
{
    typedef vector<string> strings;

    string tabs(initializer_list<const char*> p_parts)
    {
        string s;
        for (auto i = p_parts.begin(); i != p_parts.end(); ++i)
        {
            if (s.size() > 0)
            {
                s.push_back('\t');
            }
            s.insert(s.size(), *i);
        }
        return s;
    }

    typedef pair<string,uint32_t> locus;
    typedef pair<string,size_t> seq_and_count;
    typedef vector<seq_and_count> seq_and_count_list;

    struct counts
    {
        size_t nA;
        size_t nC;
        size_t nG;
        size_t nT;
        map<string,size_t> indels;

        counts()
            : nA(0), nC(0), nG(0), nT(0)
        { }

        counts(size_t p_nA, size_t p_nC, size_t p_nG, size_t p_nT, const seq_and_count_list& p_indels)
            : nA(p_nA), nC(p_nC), nG(p_nG), nT(p_nT)
        {
            for (auto itr = p_indels.begin(); itr != p_indels.end(); ++itr)
            {
                indels[itr->first] = itr->second;
            }
        }

        counts(const subtext& p_nA, const subtext& p_nC, const subtext& p_nG, const subtext& p_nT,
               const subtext& p_other)
            : nA(lexical_cast<size_t>(make_iterator_range(p_nA.first, p_nA.second))),
              nC(lexical_cast<size_t>(make_iterator_range(p_nC.first, p_nC.second))),
              nG(lexical_cast<size_t>(make_iterator_range(p_nG.first, p_nG.second))),
              nT(lexical_cast<size_t>(make_iterator_range(p_nT.first, p_nT.second)))
        {
            if (p_other == ".")
            {
                return;
            }
            vector<subtext> indel_parts;
            vector<subtext> kv;
            p_other.split(';', indel_parts);
            for (size_t i = 0; i < indel_parts.size(); ++i)
            {
                indel_parts[i].split('=', kv);
                string seq = kv[0];
                size_t cnt = lexical_cast<size_t>(make_iterator_range(kv[1].first, kv[1].second));
                indels[seq] = cnt;
            }
        }

        counts& operator+=(const counts& p_other)
        {
            nA += p_other.nA;
            nC += p_other.nC;
            nG += p_other.nG;
            nT += p_other.nT;
            for (auto itr = p_other.indels.begin(); itr != p_other.indels.end(); ++itr)
            {
                indels[itr->first] += itr->second;
            }
            return *this;
        }

        string stringify_indels() const
        {
            if (indels.size() == 0)
            {
                return ".";
            }
            else
            {
                std::string res;
                for (auto itr = indels.begin(); itr != indels.end(); ++itr)
                {
                    const std::string& seq = itr->first;
                    std::string cnt = boost::lexical_cast<std::string>(itr->second);
                    if (itr != indels.begin())
                    {
                        res.push_back(';');
                    }
                    res.insert(res.end(), seq.begin(), seq.end());
                    res.push_back('=');
                    res.insert(res.end(), cnt.begin(), cnt.end());
                }
                return res;
            }
        }
    };

    class merge_counts_command : public varoom::command
    {
    public:
        merge_counts_command(const strings& p_input_filenames,
                            const string& p_output_filename)
            : m_input_filenames(p_input_filenames),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            map<locus,counts> res;

            string chr;
            uint32_t pos;
            std::vector<subtext> other_parts;
            std::vector<subtext> other_key_val;
            string seq;

            vector<string> types{"str", "uint", "uint", "uint", "uint", "uint", "[str->uint]"};

            for (size_t n = 0; n < m_input_filenames.size(); ++n)
            {
                input_file_holder_ptr inp = files::in(m_input_filenames[n]);
                for (auto t = typed_tsv_reader(**inp, types); t.more(); ++t)
                {
                    const typed_tsv_row& r = *t;
                    chr = any_cast<string>(r[0]);
                    pos = any_cast<uint64_t>(r[1]);
                    locus loc(chr, pos);
                    res[loc] += counts(any_cast<uint64_t>(r[2]), any_cast<uint64_t>(r[3]),
                                       any_cast<uint64_t>(r[4]), any_cast<uint64_t>(r[5]),
                                       any_cast<const seq_and_count_list&>(r[6]));
                }
            }

            output_file_holder_ptr outp = files::out(m_output_filename);
            ostream& out = **outp;
            out << tabs({"chr", "pos", "nA", "nC", "nG", "nT", "indel"}) << endl;
            for (auto i = res.begin(); i != res.end(); ++i)
            {
                const string& chr = i->first.first;
                const uint32_t& pos = i->first.second;
                out << chr << '\t' << pos
                    << '\t' << i->second.nA
                    << '\t' << i->second.nC
                    << '\t' << i->second.nG
                    << '\t' << i->second.nT
                    << '\t' << i->second.stringify_indels()
                    << endl;
            }
        }

    private:
        const strings m_input_filenames;
        const string m_output_filename;
    };

    class merge_counts_factory : public command_factory
    {
    public:
        merge_counts_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("merge counts");
            opts.add_options()
                ("help,h", "produce help message")
                ("inputs", po::value<strings>(), "input filename, defaults to '-' for stdin")
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
                strings ss{"-"};
                params["inputs"] = ss;
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
            strings input_fns = p_params["inputs"];
            string output_fn = p_params["output"];
            return command_ptr(new merge_counts_command(input_fns, output_fn));
        }
    };

    bool reg = command_factory::add("merge", command_factory_ptr(new merge_counts_factory));
}
// namespace anonymous
