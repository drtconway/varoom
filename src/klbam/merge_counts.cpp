#include "varoom/command.hpp"
#include "varoom/seq/locus_stream.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/typed_tsv.hpp"

#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::seq;

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
            map<locus_id,counts> res;

            vector<string> types{"str", "uint", "uint", "uint", "uint", "uint", "[str=uint]"};

            for (size_t n = 0; n < m_input_filenames.size(); ++n)
            {
                input_file_holder_ptr inp = files::in(m_input_filenames[n]);
                locus_id prev;
                bool first = true;
                for (auto t = typed_tsv_reader(**inp, types); t.more(); ++t)
                {
                    const typed_tsv_row& r = *t;
                    locus_id loc(any_cast<const string&>(r[0]), any_cast<uint64_t>(r[1]));
                    if (first)
                    {
                        first = false;
                    }
                    else
                    {
                        if (!(prev < loc))
                        {
                            std::cerr << prev.chr() << ':' << prev.pos() << " !< " << loc.chr() << ":" << loc.pos() << std::endl;
                            throw std::runtime_error("loci out of order!");
                        }
                    }
                    prev = loc;
                    res[loc] += counts(any_cast<uint64_t>(r[2]), any_cast<uint64_t>(r[3]),
                                       any_cast<uint64_t>(r[4]), any_cast<uint64_t>(r[5]),
                                       any_cast<const seq_and_count_list&>(r[6]));
                }
            }

            vector<string> type_names{"str", "uint", "uint", "uint", "uint", "uint", "[str=uint]"};

            output_file_holder_ptr outp = files::out(m_output_filename);
            ostream& out = **outp;
            out << tabs({"chr", "pos", "nA", "nC", "nG", "nT", "indels"}) << endl;
            const tsv_column_type& ind = *tsv_column_type::get("[str=uint]");
            seq_and_count_list scl;
            string scl_str;
            for (auto i = res.begin(); i != res.end(); ++i)
            {
                const string& chr = i->first.chr();
                const uint32_t& pos = i->first.pos();

                scl.clear();
                for (auto j = i->second.indels.begin(); j != i->second.indels.end(); ++j)
                {
                    scl.push_back(*j);
                }
                ind.unmake(tsv_column_value(scl), scl_str);

                out << chr << '\t' << pos
                    << '\t' << i->second.nA
                    << '\t' << i->second.nC
                    << '\t' << i->second.nG
                    << '\t' << i->second.nT
                    << '\t' << scl_str
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
