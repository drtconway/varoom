#include "varoom/command.hpp"
#include "varoom/sam.hpp"
#include "varoom/sam/pileup.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/typed_tsv.hpp"

#include <initializer_list>

using namespace std;
using namespace varoom;

namespace // anonymous
{
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

    class pileup_holder : public varoom::sam_pileup
    {
    public:
        pileup_holder(ostream& p_out)
            : chr("<unknown>"), pos(0), out(p_out), ind(*tsv_column_type::get("[str=uint]"))
        {
        }

        string chr;
        uint32_t pos;
        ostream& out;
        const tsv_column_type& ind;
        map<string,size_t> counts;
        string other;

        virtual void output_pileup(const string& p_chr, const uint32_t& p_pos, const string& p_base, const size_t& p_count)
        {
            if (p_chr == chr && p_pos == pos)
            {
                counts[p_base] = p_count;
                return;
            }

            if (pos == 0)
            {
                chr = p_chr;
                pos = p_pos;
                counts.clear();
                counts[p_base] = p_count;
                return;
            }

            size_t nA = 0;
            size_t nC = 0;
            size_t nG = 0;
            size_t nT = 0;
            seq_and_count_list nOther;
            for (auto k = counts.begin(); k != counts.end(); ++k)
            {
                const string& seq = k->first;
                const size_t& cnt = k->second;
                if (seq == "A")
                {
                    nA = cnt;
                    continue;
                }
                if (seq == "C")
                {
                    nC = cnt;
                    continue;
                }
                if (seq == "G")
                {
                    nG = cnt;
                    continue;
                }
                if (seq == "T")
                {
                    nT = cnt;
                    continue;
                }
                nOther.push_back(pair<string,size_t>(seq, cnt));
            }

            ind.unmake(tsv_column_value(nOther), other);

            out << chr
                << '\t' << pos
                << '\t' << nA
                << '\t' << nC
                << '\t' << nG
                << '\t' << nT
                << '\t' << other
                << endl;

            chr = p_chr;
            pos = p_pos;
            counts.clear();
            counts[p_base] = p_count;
        }
    };

    class count_bases_command : public varoom::command
    {
    public:
        count_bases_command(const string& p_input_filename,
                            const string& p_output_filename)
            : m_input_filename(p_input_filename),
              m_output_filename(p_output_filename)
        {
        }

        virtual void operator()()
        {
            output_file_holder_ptr outp = files::out(m_output_filename);
            ostream& out = **outp;
            out << tabs({"chr", "pos", "nA", "nC", "nG", "nT", "indel"}) << endl;

            pileup_holder ph(out);
            input_file_holder_ptr inp = files::in(m_input_filename);
            for (sam_reader r(**inp); r.more(); ++r)
            {
                const sam_alignment& aln = *r;
                ph.add_alignment(aln.chr, aln.pos, aln.seq, aln.cigar, sam_flags::is_reverse(aln.flags));
            }
        }

    private:
        const string m_input_filename;
        const string m_output_filename;
    };

    class count_bases_factory : public command_factory
    {
    public:
        count_bases_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const std::vector<std::string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("count bases");
            opts.add_options()
                ("help,h", "produce help message")
                ("regions,r", po::value<string>(), "only include counts from regions in the named BED file")
                ("input,i", po::value<string>()->default_value("-"), "input filename, defaults to '-' for stdin")
                ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
                ;

            po::positional_options_description pos;
            pos.add("input", 1);

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

            params["input"] = vm["input"].as<string>();
            params["output"] = vm["output"].as<string>();
            if (vm.count("regions"))
            {
                params["regions"] = vm["regions"].as<string>();
            }

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
            return command_ptr(new count_bases_command(input_fn, output_fn));
        }
    };
    
    bool reg = command_factory::add("count", command_factory_ptr(new count_bases_factory));
}
// namespace anonymous

