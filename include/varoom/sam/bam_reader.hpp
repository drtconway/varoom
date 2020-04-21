#ifndef VAROOM_SAM_BAM_READER_HPP
#define VAROOM_SAM_BAM_READER_HPP

#ifndef VAROOM_SAM_SAM_ALIGNMENT_HPP
#include "varoom/sam/sam_alignment.hpp"
#endif

#include <boost/lexical_cast.hpp>

#define BOOST_LOG_DYN_LINK 1
#include <boost/log/trivial.hpp>

namespace htslib
{
#include <htslib/sam.h>
}
// namespace htslib

namespace varoom
{
    class bam_reader
    {
    public:
        bam_reader(const std::string& p_filename, bool p_verbose = false)
            : m_filename(p_filename), m_fp(NULL), m_hdr(NULL), m_aln(NULL), m_more(true), m_verbose(p_verbose), m_count(0)
        {
            m_fp = htslib::hts_open(m_filename.c_str(), "r");
            m_hdr = htslib::sam_hdr_read(m_fp);
            m_aln = htslib::bam_init1();

            next();
        }

        ~bam_reader()
        {
            if (m_aln)
            {
                htslib::bam_destroy1(m_aln);
            }
            if (m_fp)
            {
                htslib::sam_close(m_fp);
            }
        }

        bool more() const
        {
            return m_more;
        }

        const sam_alignment& operator*() const
        {
            return m_cur;
        }

        void operator++()
        {
            static const size_t M = (1 << 20) - 1;
            next();
            ++m_count;
            if (m_verbose && more() && (m_count & M) == 0)
            {
                BOOST_LOG_TRIVIAL(info) << m_filename << ": " << m_count;
            }
        }

    private:
        void next()
        {
            int res = htslib::sam_read1(m_fp, m_hdr, m_aln);
            if (res <= 0)
            {
                m_more = false;
                return;
            }

            uint32_t len = m_aln->core.l_qseq;
            const uint8_t* q = bam_get_seq(m_aln);

            m_cur.flags = m_aln->core.flag;
            if (m_aln->core.tid >= 0)
            {
                m_cur.chr = m_hdr->target_name[m_aln->core.tid];
            }
            else
            {
                m_cur.chr.clear();
                m_cur.chr.push_back('*');
            }
            m_cur.pos = m_aln->core.pos + 1;
            m_cur.mapq = m_aln->core.qual;

            m_cur.cigar.clear();
            const auto cigar = bam_get_cigar(m_aln);
            std::string n_str;
            for (size_t i = 0; i < m_aln->core.n_cigar; i++)
            {
                const int op = bam_cigar_op(cigar[i]);
                const char op_ch = BAM_CIGAR_STR[op];
                const uint32_t ol = bam_cigar_oplen(cigar[i]);

                n_str = boost::lexical_cast<std::string>(ol);
                m_cur.cigar.insert(m_cur.cigar.end(), n_str.begin(), n_str.end());
                m_cur.cigar.push_back(op_ch);
            }

            if (m_aln->core.mtid >= 0)
            {
                m_cur.mate_chr = m_hdr->target_name[m_aln->core.mtid];
            }
            else
            {
                m_cur.mate_chr.clear();
                m_cur.mate_chr.push_back('*');
            }
            m_cur.mate_pos = m_aln->core.mpos + 1;
            m_cur.tlen = m_aln->core.isize;

            m_cur.seq.clear();
            m_cur.qual.clear();
            for (size_t i = 0; i < len; ++i)
            {
                m_cur.seq.push_back(htslib::seq_nt16_str[bam_seqi(q, i)]);
                m_cur.qual.push_back(q[i]);
            }
        }

        std::string m_filename;
        htslib::samFile* m_fp;
        htslib::bam_hdr_t* m_hdr;
        htslib::bam1_t* m_aln;
        bool m_more;
        bool m_verbose;
        size_t m_count;
        sam_alignment m_cur;
    };
}
// namespace varoom

#endif // VAROOM_SAM_BAM_READER_HPP
