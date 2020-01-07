#ifndef VAROOM_UTIL_RANGES_HPP
#define VAROOM_UTIL_RANGES_HPP

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef VAROOM_UTIL_RANK_SET_HPP
#include "varoom/util/rank_set.hpp"
#endif

namespace varoom
{
    class ranges
    {
    public:
        using position = size_t;
        using range = std::pair<position,position>;
        using range_id = size_t;
        using internal_range_id = size_t;

        ranges()
            : m_next_range_id(0)
        {
        }

        range_id insert(const range& p_range)
        {
            range_id n = m_next_range_id++;
            std::vector<range> overlaps;
            overlapping_ranges(p_range.first, p_range.second, overlaps);
            range r = p_range;

            //std::cerr << "overlaps.size() = " << overlaps.size() << std::endl;

            for (size_t i = 0; i < overlaps.size(); ++i)
            {
                range& o = overlaps[i];
                //std::cerr << i
                //          << '\t' << o.first
                //          << '\t' << o.second
                //          << std::endl;
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
            return n;
        }

        void ranges_at(position p_begin, position p_end, std::vector<range_id>& p_ranges) const
        {
            p_ranges.clear();
            size_t r0 = std::min(m_begins.rank(p_begin), m_ends.rank(p_begin+1));
            size_t r1 = std::max(m_begins.rank(p_end), m_ends.rank(p_end));
            for (size_t i = r0; i < r1; ++i)
            {
                position p0 = m_begins.select(i);
                auto itr = m_index.find(p0);
                assert(itr != m_index.end());
                p_ranges.insert(p_ranges.end(), itr->second.begin(), itr->second.end());
            }
            std::sort(p_ranges.begin(), p_ranges.end());
            p_ranges.erase(std::unique(p_ranges.begin(), p_ranges.end()), p_ranges.end());
        }

        void ranges_at(position p_pos, std::vector<range_id>& p_ranges) const
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
    };
}
// namespace varoom

#endif // VAROOM_UTIL_RANGES_HPP
