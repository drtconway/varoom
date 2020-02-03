#ifndef VAROOM_UTIL_RANGES_HPP
#define VAROOM_UTIL_RANGES_HPP

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>
#include <sdsl/vlc_vector.hpp>

#ifndef VAROOM_UTIL_RANK_SET_HPP
#include "varoom/util/rank_set.hpp"
#endif

namespace varoom
{
    namespace detail
    {
    }
    // namespace detail

    class ranges_builder
    {
    public:
        using position = size_t;
        using range = std::pair<position,position>;
        using range_id = uint32_t;

        ranges_builder()
            : m_next_range_id(0)
        {
        }

        size_t size() const
        {
            return m_next_range_id;
        }

        ranges_builder& insert(const range& p_range)
        {
            range_id n = m_next_range_id++;
            m_toc[n] = p_range;

            std::vector<range> overlaps;
            overlapping_ranges(p_range.first, p_range.second, overlaps);
            range r = p_range;
            for (size_t i = 0; i < overlaps.size(); ++i)
            {
                range& o = overlaps[i];
                if (r.first < o.first)
                {
                    m_begins.insert(r.first);
                    m_ends.insert(o.first);
                    m_index[r.first].push_back(n);
                    r.first = o.first;
                }
                if (r.first > o.first)
                {
                    // split o to separate the part before r.
                    m_ends.insert(r.first);
                    m_begins.insert(r.first);
                    m_index[r.first] = m_index[o.first];
                    o.first = r.first;
                }
                assert(r.first == o.first);
                if (r.second < o.second)
                {
                    // split o to separate the part after r.
                    m_ends.insert(r.second);
                    m_begins.insert(r.second);
                    m_index[r.second] = m_index[o.first];
                    r.first = r.second;
                }
                else
                {
                    r.first = o.second;
                }
                m_index[o.first].push_back(n);
            }
            if (r.first < r.second)
            {
                m_begins.insert(r.first);
                m_ends.insert(r.second);
                m_index[r.first].push_back(n);
            }

            return *this;
        }

        void ranges_at(position p_begin, position p_end, std::vector<range>& p_ranges) const
        {
            p_ranges.clear();
            size_t r0 = std::min(m_begins.rank(p_begin), m_ends.rank(p_begin+1));
            size_t r1 = std::max(m_begins.rank(p_end), m_ends.rank(p_end));
            for (size_t i = r0; i < r1; ++i)
            {
                position p0 = m_begins.select(i);
                auto itr = m_index.find(p0);
                assert(itr != m_index.end());
                for (auto jtr = itr->second.begin(); jtr != itr->second.end(); ++jtr)
                {
                    p_ranges.push_back(m_toc.find(*jtr)->second);
                }
            }
            std::sort(p_ranges.begin(), p_ranges.end());
            p_ranges.erase(std::unique(p_ranges.begin(), p_ranges.end()), p_ranges.end());
        }

        void ranges_at(position p_pos, std::vector<range>& p_ranges) const
        {
            ranges_at(p_pos, p_pos+1, p_ranges);
        }

        const std::vector<range_id>& ranges_at_segment(position p_pos) const
        {
            auto itr = m_index.find(p_pos);
            if (itr == m_index.end())
            {
                throw std::runtime_error("no segment at position");
            }
            return itr->second;
        }

        void overlapping_ranges(position p_begin, position p_end, std::vector<range>& p_res) const
        {
            size_t r0 = std::min(m_begins.rank(p_begin), m_ends.rank(p_begin+1));
            size_t r1 = std::max(m_begins.rank(p_end), m_ends.rank(p_end));
            for (size_t i = r0; i < r1; ++i)
            {
                range x(m_begins.select(i), m_ends.select(i));
                p_res.push_back(x);
            }
        }

        const varoom::rank_set<position>& begins() const
        {
            return m_begins;
        }

        const varoom::rank_set<position>& ends() const
        {
            return m_ends;
        }

        const std::unordered_map<position,std::vector<range_id>>& index() const
        {
            return m_index;
        }

    private:
        bool in_range(position p_pos) const
        {
            size_t r0 = m_begins.rank(p_pos);
            size_t r1 = m_ends.rank(p_pos);
            return r0 != r1;
        }

        range_id m_next_range_id;
        varoom::rank_set<position> m_begins;
        varoom::rank_set<position> m_ends;
        std::unordered_map<position,std::vector<range_id>> m_index;
        std::unordered_map<range_id,range> m_toc;
    };

    class ranges
    {
    public:
        using position = size_t;
        using range = std::pair<position,position>;
        using range_id = uint32_t;
        using interval = std::pair<uint32_t,uint32_t>;

        ranges()
        {
        }

        void ranges_at(position p_begin, position p_end, std::vector<range>& p_ranges) const
        {
            p_ranges.clear();
            size_t r0 = std::min(begins_rank(p_begin), ends_rank(p_begin+1));
            size_t r1 = std::max(begins_rank(p_end), ends_rank(p_end));
            for (size_t r = r0; r < r1; ++r)
            {
                size_t i0 = m_toc_starts_select(r+1);
                size_t i1 = m_toc_starts_select(r+2);
                for (size_t i = i0; i < i1;)
                {
                    size_t br = r - m_toc_entries[i++];
                    size_t er = r + m_toc_entries[i++];
                    position b = m_begins_select(br + 1);
                    position e = m_ends_select(er + 1);
                    p_ranges.push_back(range(b, e));
                }
            }
            std::sort(p_ranges.begin(), p_ranges.end());
            p_ranges.erase(std::unique(p_ranges.begin(), p_ranges.end()), p_ranges.end());
        }

        void ranges_at(position p_pos, std::vector<range>& p_ranges) const
        {
            ranges_at(p_pos, p_pos+1, p_ranges);
        }

        void make(const ranges_builder& p_builder)
        {
            make(p_builder.begins(), m_begins, m_begins_rank, m_begins_select);
            make(p_builder.ends(), m_ends, m_ends_rank, m_ends_select);

            const std::unordered_map<size_t,std::vector<uint32_t>>& idx = p_builder.index();

            // Traverse the index (a mapping from segment start position -> list of range_id)
            // and compile a mapping from range_id to list of segment rank.
            //
            std::unordered_map<range_id,std::vector<size_t>> w;
            for (auto itr = idx.begin(); itr != idx.end(); ++itr)
            {
                size_t r = begins_rank(itr->first);
                const std::vector<uint32_t>& v = itr->second;
                for (size_t i = 0; i < v.size(); ++i)
                {
                    w[v[i]].push_back(r);
                }
            }

            std::map<size_t,std::vector<std::pair<uint32_t,uint32_t>>> jdx;
            size_t toc_count = 0;
            for (auto itr = w.begin(); itr != w.end(); ++itr)
            {
                std::vector<size_t> v = itr->second;
                std::sort(v.begin(), v.end());
                uint32_t db = 0;
                uint32_t de = v.size() - 1;
                for (size_t i = 0; i < v.size(); ++i, ++db, --de)
                {
                    jdx[v[i]].push_back(std::make_pair(db, de));
                    toc_count += 2;
                }
            }

            sdsl::bit_vector toc_starts(toc_count + 1, 0);
            sdsl::int_vector<32> toc_entries(toc_count, 0);
            size_t j = 0;
            for (auto itr = jdx.begin(); itr != jdx.end(); ++itr)
            {
                toc_starts[j] = 1;

                std::vector<std::pair<uint32_t,uint32_t>> v = itr->second;
                std::sort(v.begin(), v.end());
                for (auto jtr = v.begin(); jtr != v.end(); ++jtr)
                {
                    toc_entries[j++] = jtr->first;
                    toc_entries[j++] = jtr->second;
                }
            }
            toc_starts[j] = 1;

            m_toc_starts.swap(toc_starts);
            m_toc_starts_select = sdsl::bit_vector::select_1_type(&m_toc_starts);

            sdsl::vlc_vector<> compressed_toc_entries(toc_entries);
            m_toc_entries.swap(compressed_toc_entries);
        }

        void save(std::ostream& p_out) const
        {
            m_begins.serialize(p_out);
            m_ends.serialize(p_out);
            m_toc_starts.serialize(p_out);
            m_toc_entries.serialize(p_out);
        }

        void load(std::istream& p_in)
        {
            m_begins.load(p_in);
            m_begins_rank = sdsl::sd_vector<>::rank_1_type(&m_begins);
            m_begins_select = sdsl::sd_vector<>::select_1_type(&m_begins);
            m_ends.load(p_in);
            m_ends_rank = sdsl::sd_vector<>::rank_1_type(&m_ends);
            m_ends_select = sdsl::sd_vector<>::select_1_type(&m_ends);
            m_toc_starts.load(p_in);
            m_toc_starts_select = sdsl::bit_vector::select_1_type(&m_toc_starts);
            m_toc_entries.load(p_in);
        }

    private:
        static void make(const varoom::rank_set<position>& p_poss,
                  sdsl::sd_vector<>& p_bits,
                  sdsl::sd_vector<>::rank_1_type& p_rank,
                  sdsl::sd_vector<>::select_1_type& p_select)
        {
            size_t n = p_poss.size();
            position m = p_poss.select(n-1);

            sdsl::bit_vector B(m+1, 0);
            for (size_t i = 0; i < n; ++i)
            {
                position p = p_poss.select(i);
                B[p] = 1;
            }

            p_bits = sdsl::sd_vector<>(B);
            p_rank = sdsl::sd_vector<>::rank_1_type(&p_bits);
            p_select = sdsl::sd_vector<>::select_1_type(&p_bits);
        }

        size_t begins_rank(position p_pos) const
        {
            if (p_pos > m_begins.size())
            {
                p_pos = m_begins.size();
            }
            return m_begins_rank(p_pos);
        }

        size_t ends_rank(position p_pos) const
        {
            if (p_pos > m_ends.size())
            {
                p_pos = m_ends.size();
            }
            return m_ends_rank(p_pos);
        }

        sdsl::sd_vector<> m_begins;
        sdsl::sd_vector<>::rank_1_type m_begins_rank;
        sdsl::sd_vector<>::select_1_type m_begins_select;
        sdsl::sd_vector<> m_ends;
        sdsl::sd_vector<>::rank_1_type m_ends_rank;
        sdsl::sd_vector<>::select_1_type m_ends_select;
        sdsl::bit_vector m_toc_starts;
        sdsl::bit_vector::select_1_type m_toc_starts_select;
        sdsl::vlc_vector<> m_toc_entries;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_RANGES_HPP
