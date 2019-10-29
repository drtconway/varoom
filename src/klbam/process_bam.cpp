#include "varoom/command.hpp"
#include "varoom/sam.hpp"
#include "varoom/sam/pileup.hpp"
#include "varoom/util/files.hpp"

#include <map>

using namespace std;
using namespace varoom;

namespace // anonymous
{
    typedef pair<uint32_t,uint32_t> range;
    typedef vector<range> ranges;
    typedef map<string,ranges> regions;

    bool contains(const range& p_lhs, const uint32_t& p_x)
    {
        return p_lhs.first <= p_x && p_x < p_lhs.second;
    }

    bool overlaps(const range& p_lhs, const range& p_rhs)
    {
        return contains(p_lhs, p_rhs.first)
            || contains(p_lhs, p_rhs.second)
            || contains(p_rhs, p_lhs.first)
            || contains(p_rhs, p_lhs.second);
    }

    class pileup_dest : public sam_pileup
    {
    public:
        pileup_dest(ostream& p_out, const regions& p_regions)
            : m_out(p_out), m_regions(p_regions)
        {
        }

        virtual void output_pileup(const std::string& p_chr, const std::uint32_t& p_pos, const std::string& p_base, const size_t& p_count)
        {
            auto i = m_regions.find(p_chr);
            if (i == m_regions.end())
            {
                return;
            }
            bool found = false;
            for (auto j = i->second.begin(); j != i->second.end(); ++j)
            {
                if (contains(*j, p_pos))
                {
                    found = true;
                    break;
                }
            }
            if (found)
            {
                m_out << p_chr << '\t' << p_pos << '\t' << p_base << '\t' << p_count << endl;
            }
        }

    private:
        ostream& m_out;
        const regions& m_regions;
    };

    class process_bam : public command
    {
    public:
        process_bam(const string& p_src_filename, const string& p_dest_filename, const regions& p_regions)
            : m_src_filename(p_src_filename), m_dest_filename(p_dest_filename), m_regions(p_regions)
        {
        }

        virtual void operator()()
        {
            input_file_holder_ptr inp = files::in(m_src_filename);
            output_file_holder_ptr outp = files::out(m_dest_filename);

            pileup_dest pile(**outp, m_regions);

            string xx;
            string prevChr;
            range prevRange;
            for (sam_reader r(**inp); r.more(); ++r)
            {
                const sam_alignment& aln = *r;

                if (aln.chr != xx)
                {
                    cerr << aln.chr << endl;
                    xx = aln.chr;
                }

                range v = range(aln.pos, aln.pos + aln.seq.size());
                if (aln.chr == prevChr && overlaps(prevRange, v))
                {
                    pile.add_alignment(aln.chr, aln.pos, aln.seq, aln.cigar, sam_flags::is_reverse(aln.flags));
                    continue;
                }

                auto j = m_regions.find(aln.chr);
                if (j == m_regions.end())
                {
                    continue;
                }
                for (auto k = j->second.begin(); k != j->second.end(); ++k)
                {
                    if (overlaps(*k, v))
                    {
                        pile.add_alignment(aln.chr, aln.pos, aln.seq, aln.cigar, sam_flags::is_reverse(aln.flags));
                        prevChr = aln.chr;
                        prevRange = *k;
                        break;
                    }
                }
            }

            pile.end();
        }

    private:
        const std::string m_src_filename;
        const std::string m_dest_filename;
        const regions m_regions;
    };
}
// namespace anonymous

int main(int argc, const char** argv)
{
    regions R;
    R["18"].push_back(range(60795804, 60796044));
    R["18"].push_back(range(60985260, 60985920));

    process_bam P(argv[1], argv[2], R);
    P();
}
