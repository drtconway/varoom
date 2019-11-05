#include "varoom/command.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/tsv.hpp"

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

    typedef map<string,size_t> base_counts;
    typedef map<uint32_t,base_counts> pos_base_counts;
    typedef map<string,pos_base_counts> chr_pos_base_counts;

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
            chr_pos_base_counts counts;
            string chr;
            uint32_t pos;
            string seq;
            size_t cnt;
            for (size_t n = 0; n < m_input_filenames.size(); ++n)
            {
                input_file_holder_ptr inp = files::in(m_input_filenames[n]);
                for (auto t = tsv_reader(**inp); t.more(); ++t)
                {
                    const tsv_row& r = *t;
                    chr = r[0];
                    pos = lexical_cast<uint32_t>(make_iterator_range(r[1].first, r[1].second));
                    seq = r[2];
                    cnt = lexical_cast<uint32_t>(make_iterator_range(r[3].first, r[3].second));
                    counts[chr][pos][seq] += cnt;
                }
            }

            output_file_holder_ptr outp = files::out(m_output_filename);
            ostream& out = **outp;
            out << tabs({"chr", "pos", "seq", "count"}) << endl;
            for (auto i = counts.begin(); i != counts.end(); ++i)
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
        const strings m_input_filenames;
        const string m_output_filename;
    };

    command_ptr merge_counts_factory(const command_parameters& p_params)
    {
        strings input_fns({"-"});
        if (p_params.count("input"))
        {
            input_fns = p_params["input"].as<strings>();
        }
        string output_fn = p_params["output"].as<string>();
        return command_ptr(new merge_counts_command(input_fns, output_fn));
    }

    command_options_ptr merge_counts_options()
    {
        namespace po = boost::program_options;

        command_options_ptr opts_ptr(new command_options("merge counts"));
        (*opts_ptr).add_options()
            ("help,h", "produce help message")
            ("input,i", po::value<strings>(), "input filenames, defaults to '-' for stdin")
            ("output,o", po::value<string>()->default_value("-"), "output filename, defaults to '-' for stdout")
            ;
        return opts_ptr;
    }

    bool reg = command_factory::add("merge", merge_counts_factory, merge_counts_options());
}
// namespace anonymous
