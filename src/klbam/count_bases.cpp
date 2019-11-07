#include "varoom/command.hpp"
#include "varoom/sam.hpp"
#include "varoom/sam/pileup.hpp"
#include "varoom/util/files.hpp"

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

    class pileup_holder : public varoom::sam_pileup
    {
    public:
        typedef map<string,size_t> base_counts;
        typedef map<uint32_t,base_counts> pos_base_counts;
        typedef map<string,pos_base_counts> chr_pos_base_counts;

        chr_pos_base_counts piles;

        virtual void output_pileup(const string& p_chr, const uint32_t& p_pos, const string& p_base, const size_t& p_count)
        {
            piles[p_chr][p_pos][p_base] = p_count;
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
            pileup_holder ph;

            input_file_holder_ptr inp = files::in(m_input_filename);
            for (sam_reader r(**inp); r.more(); ++r)
            {
                const sam_alignment& aln = *r;
                ph.add_alignment(aln.chr, aln.pos, aln.seq, aln.cigar, sam_flags::is_reverse(aln.flags));
            }

            output_file_holder_ptr outp = files::out(m_output_filename);
            ostream& out = **outp;
            out << tabs({"chr", "pos", "seq", "count"}) << endl;
            for (auto i = ph.piles.begin(); i != ph.piles.end(); ++i)
            {
                const string& chr = i->first;
                for (auto j = i->second.begin(); j != i->second.end(); ++j)
                {
                    const uint32_t& pos = j->first;
                    for (auto k = j->second.begin(); k != j->second.end(); ++k)
                    {
                        const string& seq = k->first;
                        const size_t& cnt = k->second;
                        out << chr
                            << '\t' << pos
                            << '\t' << seq
                            << '\t' << cnt
                            << endl;
                    }
                }
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

