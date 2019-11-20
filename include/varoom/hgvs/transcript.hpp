#ifndef VAROOM_HGVS_TRANSCRIPT_HPP
#define VAROOM_HGVS_TRANSCRIPT_HPP

#ifndef VAROOM_HGVS_HGVS_POSITIONS_HPP
#include "varoom/hgvs/hgvs_positions.hpp"
#endif

#ifndef VAROOM_HGVS_HGVSC_LOCUS_HPP
#include "varoom/hgvs/hgvsc_locus.hpp"
#endif

#include <algorithm>
#include <iostream>
#include <vector>

namespace varoom
{
    namespace hgvs
    {
        enum tx_strand { POS, NEG };

        typedef std::pair<genomic_locus,genomic_locus> genomic_exon;
        typedef std::pair<tx_position,tx_position> tx_exon;

        class transcript
        {
        public:
            transcript(const std::string& p_accession,
                       const std::string& p_chr,
                       const tx_strand& p_strand,
                       const genomic_locus& p_tx_begin, const genomic_locus& p_tx_end,
                       const genomic_locus& p_cds_begin, const genomic_locus& p_cds_end,
                       const std::vector<genomic_exon>& p_exons)
                : m_accession(p_accession), m_chr(p_chr), m_strand(p_strand),
                  m_tx_begin(p_tx_begin), m_tx_end(p_tx_end),
                  m_cds_begin(p_cds_begin), m_cds_end(p_cds_end),
                  m_g_exons(p_exons)
            {
                  make_t_exons();
            }

            hgvsc_locus to_hgvsc_locus(const hgvsg_locus& p_gpos) const
            {
                genomic_locus g0(static_cast<uint64_t>(p_gpos) - 1);

                if (g0 < m_tx_begin || g0 >= m_tx_end)
                {
                    return intergenic(g0);
                }
                if (g0 < m_cds_begin)
                {
                    size_t exNo = exon_of(g0);
                    return upstreamUtr(g0, exNo);
                }
                if (g0 >= m_cds_end)
                {
                    size_t exNo = exon_of(g0);
                    return downstreamUtr(g0, exNo);
                }
                size_t exNo = exon_of(g0);
                if (contains(m_g_exons[exNo], g0))
                {
                    throw std::runtime_error("coding not implemented");
                    //return coding(g0, exNo);
                }
                throw std::runtime_error("intronic not implemented");
                //return intronic(g0, exNo);

            }

            hgvsc_locus intergenic(genomic_locus& p_gpos) const
            {
                switch (m_strand)
                {
                    case POS:
                    {
                        if (p_gpos < m_tx_begin)
                        {
                            int64_t d = static_cast<uint64_t>(m_tx_begin - p_gpos);
                            tx_position p0 = m_t_exons.front().first - tx_position(d);
                            hgvsc_position p((static_cast<int64_t>(p0) - 1));
                            return hgvsc_locus(UTR5, p, hgvsc_relative_position(0));
                        }
                        if (p_gpos >= m_tx_end)
                        {
                            int64_t d = static_cast<uint64_t>(p_gpos - m_tx_end);
                            tx_position p0 = m_t_exons.back().second - m_t_cds_end + tx_position(d);
                            // Because intervals are half-open, we don't need to shift the position.
                            hgvsc_position p(static_cast<int64_t>(p0));
                            return hgvsc_locus(UTR3, p, hgvsc_relative_position(0));
                        }
                        throw std::runtime_error("internal logic error");
                    }

                    case NEG:
                    {
                        if (p_gpos < m_tx_begin)
                        {
                            int64_t d = static_cast<uint64_t>(m_tx_begin - p_gpos);
                            tx_position p0 = (m_t_exons.back().second - m_t_cds_end) + tx_position(d);
                            hgvsc_position p((static_cast<int64_t>(p0) + 1));
                            return hgvsc_locus(UTR3, p, hgvsc_relative_position(0));
                        }
                        if (p_gpos >= m_tx_end)
                        {
                            int64_t d = static_cast<uint64_t>(p_gpos - m_tx_end);
                            tx_position p0 = m_t_exons.front().first - tx_position(d);
                            // Because intervals are half-open, we don't need to shift the position.
                            hgvsc_position p(static_cast<int64_t>(p0));
                            return hgvsc_locus(UTR5, p, hgvsc_relative_position(0));
                        }
                        throw std::runtime_error("internal logic error");
                    }
                }
                throw std::runtime_error("internal logic error");
            }

            hgvsc_locus upstreamUtr(genomic_locus& p_gpos, const size_t p_g_exNo) const
            {
                switch (m_strand)
                {
                    case POS:
                    {
                        size_t t_exNo = gexToTex(p_g_exNo);
                        const genomic_exon& g_ex = m_g_exons[p_g_exNo];
                        const tx_exon& t_ex = m_t_exons[t_exNo];
                        if (contains(m_g_exons[p_g_exNo], p_gpos))
                        {
                            int64_t p0 = static_cast<uint64_t>(p_gpos - g_ex.first);
                            tx_position p1 = t_ex.first + tx_position(p0);
                            hgvsc_position p(static_cast<int64_t>(p1) - 1);
                            return hgvsc_locus(UTR5, p, hgvsc_relative_position(0));
                        }

                        // Intronic
                        if (p_gpos < g_ex.first)
                        {
                            uint64_t d3 = static_cast<uint64_t>(g_ex.first - p_gpos);
                            hgvsc_position p(static_cast<int64_t>(t_ex.first));
                            hgvsc_relative_position r(-d3 - 1);
                            return hgvsc_locus(UTR5, p, r);
                        }
                        if (p_gpos >= g_ex.second)
                        {
                            uint64_t d5 = static_cast<uint64_t>(p_gpos - g_ex.second);
                            hgvsc_position p(static_cast<int64_t>(t_ex.second) - 1);
                            hgvsc_relative_position r(d5);
                            return hgvsc_locus(UTR5, p, r);
                        }
                        throw std::runtime_error("internal logic error");
                    }

                    case NEG:
                    {
                        size_t t_exNo = gexToTex(p_g_exNo);
                        const genomic_exon& g_ex = m_g_exons[p_g_exNo];
                        const tx_exon& t_ex = m_t_exons[t_exNo];
                        if (contains(g_ex, p_gpos))
                        {
                            int64_t p0 = static_cast<uint64_t>(g_ex.second - p_gpos);
                            tx_position p1 = t_ex.first + tx_position(p0) - m_t_cds_end;
                            hgvsc_position p(static_cast<int64_t>(p1) + 1);
                            return hgvsc_locus(UTR3, p, hgvsc_relative_position(0));
                        }

                        // Intronic
                        if (p_gpos < g_ex.first)
                        {
                            uint64_t d3 = static_cast<uint64_t>(g_ex.first - p_gpos);
                            hgvsc_position p(static_cast<int64_t>(t_ex.second - m_t_cds_end));
                            hgvsc_relative_position r(d3 + 1);
                            return hgvsc_locus(UTR3, p, r);
                        }
                        if (p_gpos >= g_ex.second)
                        {
                            uint64_t d5 = static_cast<uint64_t>(p_gpos - g_ex.second);
                            hgvsc_position p(static_cast<int64_t>(t_ex.first - m_t_cds_end) + 1);
                            hgvsc_relative_position r(-d5);
                            return hgvsc_locus(UTR3, p, r);
                        }
                        throw std::runtime_error("internal logic error");
                    }
                }
                throw std::runtime_error("internal logic error");
            }

            hgvsc_locus downstreamUtr(genomic_locus& p_gpos, const size_t p_g_exNo) const
            {
                switch (m_strand)
                {
                    case POS:
                    {
                        size_t t_exNo = gexToTex(p_g_exNo);
                        const genomic_exon& g_ex = m_g_exons[p_g_exNo];
                        const tx_exon& t_ex = m_t_exons[t_exNo];
                        if (contains(m_g_exons[p_g_exNo], p_gpos))
                        {
                            int64_t p0 = static_cast<uint64_t>(p_gpos - g_ex.first);
                            tx_position p1 = t_ex.first - m_t_cds_end + tx_position(p0);
                            hgvsc_position p(static_cast<int64_t>(p1));
                            return hgvsc_locus(UTR3, p, hgvsc_relative_position(0));
                        }

                        // Intronic
                        if (p_gpos < g_ex.first)
                        {
                            uint64_t d3 = static_cast<uint64_t>(g_ex.first - p_gpos);
                            hgvsc_position p(static_cast<int64_t>(t_ex.first - m_t_cds_end) + 1);
                            hgvsc_relative_position r(-d3 - 1);
                            return hgvsc_locus(UTR3, p, r);
                        }
                        if (p_gpos >= g_ex.second)
                        {
                            uint64_t d5 = static_cast<uint64_t>(p_gpos - g_ex.second);
                            hgvsc_position p(static_cast<int64_t>(t_ex.second - m_t_cds_end));
                            hgvsc_relative_position r(d5);
                            return hgvsc_locus(UTR3, p, r);
                        }
                        throw std::runtime_error("internal logic error");
                    }

                    case NEG:
                    {
                        size_t t_exNo = gexToTex(p_g_exNo);
                        const genomic_exon& g_ex = m_g_exons[p_g_exNo];
                        const tx_exon& t_ex = m_t_exons[t_exNo];
                        if (contains(g_ex, p_gpos))
                        {
                            int64_t p0 = static_cast<uint64_t>(g_ex.second - p_gpos);
                            tx_position p1 = t_ex.first + tx_position(p0);
                            hgvsc_position p(static_cast<int64_t>(p1));
                            return hgvsc_locus(UTR5, p, hgvsc_relative_position(0));
                        }

                        // Intronic
                        if (p_gpos < g_ex.first)
                        {
                            uint64_t d3 = static_cast<uint64_t>(g_ex.first - p_gpos);
                            hgvsc_position p(static_cast<int64_t>(t_ex.second) - 1);
                            hgvsc_relative_position r(d3 + 1);
                            return hgvsc_locus(UTR5, p, r);
                        }
                        if (p_gpos >= g_ex.second)
                        {
                            uint64_t d5 = static_cast<uint64_t>(p_gpos - g_ex.second);
                            hgvsc_position p(static_cast<int64_t>(t_ex.first));
                            hgvsc_relative_position r(-d5);
                            return hgvsc_locus(UTR5, p, r);
                        }
                        throw std::runtime_error("internal logic error");
                    }
                }
                throw std::runtime_error("internal logic error");
            }

            size_t exon_of(const genomic_locus& p_gpos) const
            {
                for (size_t i = 0; i < m_g_exons.size(); ++i)
                {
                    if (p_gpos < m_g_exons[i].first)
                    {
                        if (i == 0)
                        {
                            return i;
                        }
                        uint64_t d5 = static_cast<uint64_t>(p_gpos - m_g_exons[i - 1].second) + 1;
                        uint64_t d3 = static_cast<uint64_t>(m_g_exons[i].first - p_gpos);
                        if (d5 < d3)
                        {
                            return i - 1;
                        }
                        else
                        {
                            return i;
                        }
                    }
                    if (m_g_exons[i].first <= p_gpos && p_gpos < m_g_exons[i].second)
                    {
                        return i;
                    }
                }
                return m_g_exons.size() - 1;
            }

            size_t gexToTex(const size_t p_g_exNo) const
            {
                if (m_strand == POS)
                {
                    return p_g_exNo;
                }
                else
                {
                    return m_g_exons.size() - p_g_exNo - 1;
                }
            }
#if 0
            hgvsc_locus to_txpos(const hgvsg_locus& p_gpos) const
            {
                genomic_locus gpos = static_cast<genomic_locus>(p_gpos);
                if (gpos < m_tx_begin || gpos >= m_tx_end)
                {
                    return intergenic(gpos);
                }
                size_t exNo = exon_of(gpos);
                if (gpos < m_cds_begin)
                {
                    return upstreamUtr(gpos, exNo);
                }
                if (gpos >= m_cds_end)
                {
                    return downstreamUtr(gpos, exNo);
                }
                if (contains(exNo, gpos))
                {
                    return coding(gpos, exNo);
                }
                return intronic(gpos, exNo);
            }

            hgvsg_locus to_gpos(const hgvsc_locus& p_txpos) const
            {
            }
#endif

        private:
#if 0
            hgvsc_locus intergenic(const genomic_locus& p_gpos, const size_t& p_exNo) const
            {
                switch (m_strand)
                {
                    case POS:
                    {
                        if (p_gpos < m_tx_begin)
                        {
                            hgvsc_locus res;
                            res.zone = UTR5;
                            res.tx_pos = m_tx_exons.front().first() + (m_tx_begin() - p_gpos());
                            res.rel_pos = 0;
                            return res;
                        }
                        else
                        {
                            hgvsc_locus res;
                            res.zone = UTR3;
                            res.tx_pos = m_tx_exons.back().second() + (p_gpos() - m_tx_end()) + 1;
                            res.rel_pos = 0;
                            return res;
                        }
                    }
                    case NEG:
                    {
                        if (p_gpos < m_tx_begin)
                        {
                            hgvsc_locus res;
                            res.zone = UTR3;
                            res.tx_pos = m_tx_exons.back().second() + (m_tx_begin() - p_gpos()) + 1;
                            res.rel_pos = 0;
                            return res;
                        }
                        else
                        {
                            hgvsc_locus res;
                            res.zone = UTR5;
                            res.tx_pos = m_tx_exons.front().first() + (m_tx_end() - p_gpos());
                            res.rel_pos = 0;
                            return res;
                        }
                    }
                }
            }

            bool contains(const size_t& p_exon_num, const genomic_locus& p_gpos) const
            {
                const std::pair<genomic_locus,genomic_locus>& exon = m_g_exons[p_exon_num];
                return exon.first <= p_gpos && p_gpos < exon.second;
            }

            size_t exon_of(const genomic_locus& p_gpos) const
            {
                std::int64_t pos = p_gpos();
                for (size_t i = 0; i < m_g_exons.size(); ++i)
                {
                    const std::pair<genomic_locus,genomic_locus>& exon = m_g_exons[i];
                    if (pos < exon.first)
                    {
                        return nearest_exon(p_gpos, i - 1, i);
                    }
                    if (exon.first <= pos && pos < exon.second)
                    {
                        return i;
                    }
                }
                return nearest_exon(p_gpos, m_g_exons.size() - 1, m_g_exons.size();
            }

            size_t nearest_exon(const genomic_locus& p_gpos, const size_t& p_lhs, const size_t& p_rhs) const
            {
                if (p_lhs == -1)
                {
                    return p_rhs;
                }
                if (p_rhs == m_g_exons.size())
                {
                    return p_lhs;
                }
                std::int64_t pos = p_gpos();
                const std::pair<genomic_locus,genomic_locus>& lhsExon = m_g_exons[p_lhs];
                const std::pair<genomic_locus,genomic_locus>& rhsExon = m_g_exons[p_rhs];
                std::int64_t lhsDist = pos - lhsExon.second() + 1;
                std::int64_t rhsDist = rhsExon.first() - pos;
                if (lhsDist < rhsDist)
                {
                    return p_lhs;
                }
                else
                {
                    return p_rhs;
                }
            }

            bool contains(const size_t& p_exon_num, const tx_position& p_gpos) const
            {
                const std::pair<tx_position,tx_position>& exon = m_tx_exons[p_exon_num];
                return exon.first <= p_gpos && p_gpos < exon.second;
            }

            size_t exon_of(const tx_position& p_txpos) const
            {
                std::int64_t pos = p_txpos();
                for (size_t i = 0; i < m_tx_exons.size(); ++i)
                {
                    const std::pair<tx_position,tx_position>& exon = m_tx_exons[i];
                    if (pos < exon.first)
                    {
                        return nearest_exon(p_txpos, i - 1, i);
                    }
                    if (exon.first <= pos && pos < exon.second)
                    {
                        return i;
                    }
                }
                return nearest_exon(p_txpos, m_tx_exons.size() - 1, m_tx_exons.size();
            }

            size_t nearest_exon(const tx_position& p_gpos, const size_t& p_lhs, const size_t& p_rhs) const
            {
                if (p_lhs == -1)
                {
                    return p_rhs;
                }
                if (p_rhs == m_tx_exons.size())
                {
                    return p_lhs;
                }
                std::int64_t pos = p_txpos();
                const std::pair<tx_position,tx_position>& lhsExon = m_tx_exons[p_lhs];
                const std::pair<tx_position,tx_position>& rhsExon = m_tx_exons[p_rhs];
                std::int64_t lhsDist = pos - lhsExon.second() + 1;
                std::int64_t rhsDist = rhsExon.first() - pos;
                if (lhsDist < rhsDist)
                {
                    return p_lhs;
                }
                else
                {
                    return p_rhs;
                }
            }
#endif

            void make_t_exons()
            {
                typedef std::pair<int64_t,int64_t> pos_pair;

                std::vector<pos_pair> xs;
                int64_t tb = 0;
                int64_t te = 0;
                int64_t t = 0;
                for (size_t i = 0; i < m_g_exons.size(); ++i)
                {
                    const genomic_exon& ex = m_g_exons[i];
                    uint64_t l = static_cast<uint64_t>(ex.second) - static_cast<uint64_t>(ex.first);

                    int64_t b = t;
                    int64_t e = t + l;
                    //std::cout << i << '\t' << b << '\t' << e << std::endl;
                    xs.push_back(pos_pair(b, e));
                    if (contains(ex, m_cds_begin))
                    {
                        tb = b + static_cast<uint64_t>(m_cds_begin - ex.first);
                    }
                    if (contains(ex, m_cds_end))
                    {
                        te = b + static_cast<uint64_t>(m_cds_end - ex.first);
                    }
                    t += l;
                }
                //std::cout << tb << '\t' << te << std::endl;

                switch (m_strand)
                {
                    case POS:
                    {
                        for (size_t i = 0; i < xs.size(); ++i)
                        {
                            xs[i].first -= tb;
                            xs[i].second -= tb;
                        }
                        m_t_cds_begin = tx_position(tb);
                        m_t_cds_end = tx_position(te) - m_t_cds_begin;
                        break;
                    }
                    case NEG:
                    {
                        for (size_t i = 0; i < xs.size(); ++i)
                        {
                            xs[i].first = te - xs[i].first;
                            xs[i].second = te - xs[i].second;
                            std::swap(xs[i].first, xs[i].second);
                        }
                        std::reverse(xs.begin(), xs.end());
                        m_t_cds_begin = tx_position(t - te);
                        m_t_cds_end = tx_position(t - tb) - m_t_cds_begin;
                        break;
                    }
                }

                for (size_t i = 0; i < xs.size(); ++i)
                {
                    //std::cout << i << '\t' << xs[i].first << '\t' << xs[i].second << std::endl;
                    m_t_exons.push_back(tx_exon(tx_position(xs[i].first), tx_position(xs[i].second)));
                }

                m_t_tx_begin = m_t_exons.front().first;
                m_t_tx_end = m_t_exons.back().second;
            }

            static bool contains(const genomic_exon& p_exon, const genomic_locus& p_loc)
            {
                return p_exon.first <= p_loc && p_loc < p_exon.second;
            }

            const std::string m_accession;
            const std::string m_chr;
            const tx_strand m_strand;
            const genomic_locus m_tx_begin;
            const genomic_locus m_tx_end;
            const genomic_locus m_cds_begin;
            const genomic_locus m_cds_end;
            const std::vector<genomic_exon> m_g_exons;
            tx_position m_t_tx_begin;
            tx_position m_t_tx_end;
            tx_position m_t_cds_begin;
            tx_position m_t_cds_end;
            std::vector<tx_exon> m_t_exons;
        };
    }
    // namespace hgvs
}
// namespace varoom

#endif // VAROOM_HGVS_TRANSCRIPT_HPP
