#include "varoom/sam/bam_reader.hpp"
#include "varoom/command.hpp"
#include "varoom/sam.hpp"
#include "varoom/sam/pileup.hpp"
#include "varoom/seq/locus_ordering.hpp"
#include "varoom/util/files.hpp"
#include "varoom/util/text.hpp"
#include "varoom/util/tsv.hpp"
#include "varoom/util/table.hpp"
#include "varoom/klbam/types.hpp"

#include <initializer_list>
#include <tuple>

using namespace std;
using namespace boost;
using namespace varoom;
using namespace varoom::seq;
using namespace varoom::klbam;


namespace // anonymous
{
    typedef pair<uint32_t,uint32_t> region;
    typedef vector<region> region_list;
    typedef map<string,region_list> regions;

    void read_bed_file(const string& p_filename, regions& p_res)
    {
        input_file_holder_ptr inp = files::in(p_filename);
        for (tsv_reader in(**inp); in.more(); ++in)
        {
            const tsv_row& r = *in;
            region v;
            v.first = r.get<uint32_t>(1);
            v.second = r.get<uint32_t>(2);
            p_res[string(r[0])].push_back(v);
        }

        // Make sure they're sorted
        //
        for (auto itr = p_res.begin(); itr != p_res.end(); ++itr)
        {
            sort(itr->second.begin(), itr->second.end());
        }
    }

    bool contains(const regions& p_reg, const string& p_chr, uint32_t p_pos)
    {
        // The empty set is interpreted as include everything!
        //
        if (p_reg.size() == 0)
        {
            return true;
        }

        auto itr = p_reg.find(p_chr);
        if (itr == p_reg.end())
        {
            return false;
        }
        const region_list& rs = itr->second;

        auto fst = rs.begin();
        auto snd = rs.end();
        size_t count = distance(fst, snd);
        while (count > 0)
        {
            auto it = fst;
            size_t step = count / 2;
            std::advance(it, step);
            if (it->second < p_pos)
            {
                fst = ++it;
                count -= step + 1;
            }
            else
            {
                count = step;
            }
        }
        if (fst != snd)
        {
            return fst->first <= p_pos && p_pos < fst->second;
        }
        return false;

        if (0)
        {
            for (size_t i = 0; i < rs.size(); ++i)
            {
                const region& r = rs[i];
                if (r.first <= p_pos && p_pos < r.second)
                {
                    return true;
                }
            }
            return false;
        }
    }


    class pileup_holder : public varoom::sam_pileup
    {
    public:
        pileup_holder(vector<counts::tuple>& p_out, const regions& p_reg)
            : chr("<unknown>"), pos(0), out(p_out), reg(p_reg)
        {
        }

    private:
        virtual void output_pileup(const string& p_chr, const uint32_t& p_pos, const string& p_base, const size_t& p_count)
        {
            if (!contains(reg, p_chr, p_pos))
            {
                return;
            }

            if (p_chr == chr && p_pos == pos)
            {
                base_counts[p_base] = p_count;
                return;
            }

            if (pos == 0)
            {
                chr = p_chr;
                pos = p_pos;
                base_counts.clear();
                base_counts[p_base] = p_count;
                return;
            }

            counts::tuple row(chr, pos, 0, 0, 0, 0, 0, 0, 0, 0, 0, "");
            map<std::string,size_t> nOther;
            for (auto k = base_counts.begin(); k != base_counts.end(); ++k)
            {
                const string& seq = k->first;
                const size_t& cnt = k->second;
                if (seq == "A")
                {
                    std::get<counts::nA>(row) += cnt;
                    continue;
                }
                if (seq == "C")
                {
                    std::get<counts::nC>(row) += cnt;
                    continue;
                }
                if (seq == "G")
                {
                    std::get<counts::nG>(row) += cnt;
                    continue;
                }
                if (seq == "T")
                {
                    std::get<counts::nT>(row) += cnt;
                    continue;
                }
                if (seq == "N")
                {
                    std::get<counts::nN>(row) += cnt;
                    continue;
                }
                if (seq == ">")
                {
                    std::get<counts::nX0>(row) += cnt;
                    continue;
                }
                if (seq == "<")
                {
                    std::get<counts::nX1>(row) += cnt;
                    continue;
                }
                if (seq == "I")
                {
                    std::get<counts::nI>(row) += cnt;
                    continue;
                }
                if (seq == "D")
                {
                    std::get<counts::nD>(row) += cnt;
                    continue;
                }
                nOther[seq] += cnt;
            }

            other.clear();
            for (auto j = nOther.begin(); j != nOther.end(); ++j)
            {
                std::string v = j->first + std::string("=") + boost::lexical_cast<std::string>(j->second);
                if (other.size() > 0)
                {
                    other.push_back(';');
                }
                other.insert(other.end(), v.begin(), v.end());
            }
            std::get<counts::oth>(row) = std::move(other);
            out.push_back(row);

            chr = p_chr;
            pos = p_pos;
            base_counts.clear();
            base_counts[p_base] = p_count;
        }

        string chr;
        uint32_t pos;
        vector<counts::tuple>& out;
        const regions& reg;
        map<string,size_t> base_counts;
        string other;
    };

    class count_bases_command : public varoom::command
    {
    public:
        count_bases_command(const string& p_input_filename,
                            const string& p_output_filename,
                            bool p_force_sam, bool p_force_bam)
            : m_input_filename(p_input_filename),
              m_output_filename(p_output_filename),
              m_force_sam(p_force_sam), m_force_bam(p_force_bam)
        {
        }

        count_bases_command(const string& p_input_filename,
                            const string& p_output_filename,
                            const string& p_regions_filename,
                            bool p_force_sam, bool p_force_bam)
            : m_input_filename(p_input_filename),
              m_output_filename(p_output_filename),
              m_force_sam(p_force_sam), m_force_bam(p_force_bam)
        {
            read_bed_file(p_regions_filename, m_regions);
        }

        virtual void operator()()
        {
            bool use_sam = true; // if there's no filename, default to SAM.
            if (m_force_sam && m_force_bam)
            {
                cerr << "usage error: both --sam|-s and --bam|-b specificed" << endl;
                return;
            }
            if (m_force_sam || m_force_bam)
            {
                use_sam = m_force_sam;
            }
            else if (ends_with(m_input_filename, ".sam") || ends_with(m_input_filename, ".sam.gz"))
            {
                use_sam = true;
            }
            else if (ends_with(m_input_filename, ".bam"))
            {
                use_sam = false;
            }
            
            counts::table rows;
            pileup_holder ph(rows, m_regions);
            if (use_sam)
            {
                input_file_holder_ptr inp = files::in(m_input_filename);
                for (sam_reader r(**inp); r.more(); ++r)
                {
                    const sam_alignment& aln = *r;
                    if (sam_flags::is_unmapped(aln.flags))
                    {
                        continue;
                    }
                    ph.add_alignment(aln.chr, aln.pos, aln.seq, aln.cigar);
                }
            }
            else
            {
                for (bam_reader r(m_input_filename); r.more(); ++r)
                {
                    const sam_alignment& aln = *r;
                    if (sam_flags::is_unmapped(aln.flags))
                    {
                        continue;
                    }
                    ph.add_alignment(aln.chr, aln.pos, aln.seq, aln.cigar);
                }
            }

            std::sort(rows.begin(), rows.end());
            output_file_holder_ptr outp = files::out(m_output_filename);
            ostream& out = **outp;
            counts::write(out, rows, counts::labels());
        }

    private:
        static bool ends_with(const string& p_str, const string& p_suffix)
        {
            auto s = p_str.rbegin();
            auto t = p_suffix.rbegin();
            while (s != p_str.rend() && t != p_suffix.rend())
            {
                if (*s != *t)
                {
                    return false;
                }
                ++s;
                ++t;
            }
            return t == p_suffix.rend();
        }

        const string m_input_filename;
        const string m_output_filename;
        const bool m_force_sam;
        const bool m_force_bam;
        regions m_regions;
    };

    class count_bases_factory : public command_factory
    {
    public:
        count_bases_factory() {}

        virtual json parse(const command_options& p_global_opts, const command_parameters& p_globals,
                           const vector<string>& p_args) const
        {
            namespace po = boost::program_options;

            command_options opts("count bases");
            opts.add_options()
                ("help,h", "produce help message")
                ("bam,b", "force the interpretation of the input as BAM")
                ("sam,s", "force the interpretation of the input as SAM")
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
            params["sam"] = (vm.count("sam") > 0);
            params["bam"] = (vm.count("bam") > 0);

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
            bool force_sam = p_params["sam"];
            bool force_bam = p_params["bam"];
            if (p_params.count("regions"))
            {
                string regions_fn = p_params["regions"];
                return command_ptr(new count_bases_command(input_fn, output_fn, regions_fn, force_sam, force_bam));
            }
            else
            {
                return command_ptr(new count_bases_command(input_fn, output_fn, force_sam, force_bam));
            }
        }
    };
    
    bool reg = command_factory::add("count", command_factory_ptr(new count_bases_factory));
}
// namespace anonymous

