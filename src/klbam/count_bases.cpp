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

    command_ptr count_bases_factory(const command_parameters& p_params)
    {
        string input_fn = p_params["input"].as<string>();
        string output_fn = p_params["output"].as<string>();
        return command_ptr(new count_bases_command(input_fn, output_fn));
    }

    command_options_ptr count_bases_options()
    {
        namespace po = boost::program_options;

        command_options_ptr opts_ptr(new command_options("count bases"));
        (*opts_ptr).add_options()
            ("help,h", "produce help message")
            ("regions,r", po::value<string>(), "only include counts from regions in the named BED file")
            ("input,i", po::value<string>()->default_value("-"), "input filename, defaults to '-' for stdin")
            ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
            ;
        return opts_ptr;
    }

    bool reg = command_factory::add("count", count_bases_factory, count_bases_options());
}
// namespace anonymous

