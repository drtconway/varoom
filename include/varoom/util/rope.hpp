#ifndef VAROOM_UTIL_ROPE_HPP
#define VAROOM_UTIL_ROPE_HPP

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace varoom
{
    class rope
    {
    public:

        using str_range = std::pair<std::string::const_iterator,std::string::const_iterator>;
        using str_ranges = std::vector<str_range>;

        class node
        {
        public:
            node(size_t p_size)
                : m_size(p_size)
            {
            }

            virtual ~node() {}

            size_t size() const
            {
                return m_size;
            }

            virtual char operator[](size_t p_idx) const = 0;

            virtual void str(size_t p_begin, size_t p_end, str_ranges& p_res) const = 0;

        private:
            const size_t m_size;
        };
        using node_ptr = std::shared_ptr<node>;

        class atom : public node
        {
        public:
            atom(const std::string& p_atom)
                : node(p_atom.size()), m_atom(p_atom)
            {
            }

            virtual char operator[](size_t p_idx) const
            {
                return m_atom[p_idx];
            }

            virtual void str(size_t p_begin, size_t p_end, str_ranges& p_res) const
            {
                p_res.push_back(str_range(m_atom.begin() + p_begin, m_atom.begin() + p_end));
            }

        private:
            const std::string m_atom;
        };

        class concat : public node
        {
        public:
            concat(node_ptr p_lhs, node_ptr p_rhs)
                : node(p_lhs->size() + p_rhs->size()), m_lhs(p_lhs), m_rhs(p_rhs)
            {
            }

            virtual char operator[](size_t p_idx) const
            {
                size_t lz = m_lhs->size();
                if (p_idx < lz)
                {
                    return (*m_lhs)[p_idx];
                }
                else
                {
                    return (*m_rhs)[p_idx - lz];
                }
            }

            virtual void str(size_t p_begin, size_t p_end, str_ranges& p_res) const
            {
                size_t lz = m_lhs->size();
                if (p_begin < lz && p_end <= lz)
                {
                    m_lhs->str(p_begin, p_end, p_res);
                    return;
                }
                if (p_begin >= lz)
                {
                    m_rhs->str(p_begin - lz, p_end - lz, p_res);
                    return;
                }
                m_lhs->str(p_begin, lz, p_res);
                m_rhs->str(0, p_end - lz, p_res);
            }

        private:
            const node_ptr m_lhs;
            const node_ptr m_rhs;
        };

        
        class substr : public node
        {
        public:
            substr(node_ptr p_parent, size_t p_begin, size_t p_end)
                : node(p_end - p_begin), m_parent(p_parent), m_begin(p_begin), m_end(p_end)
            {
            }

            virtual char operator[](size_t p_idx) const
            {
                return (*m_parent)[m_begin + p_idx];
            }

            virtual void str(size_t p_begin, size_t p_end, str_ranges& p_res) const
            {
                m_parent->str(m_begin + p_begin, m_begin + p_end, p_res);
            }

        private:
            const node_ptr m_parent;
            const size_t m_begin;
            const size_t m_end;
        };

        class block : public node
        {
        public:
            block(size_t p_size, size_t p_block_bits, std::function<std::string(size_t)> p_func)
                : node(p_size), m_block_bits(p_block_bits), m_func(p_func)
            {
            }

            virtual char operator[](size_t p_idx) const
            {
                size_t b = block_num(p_idx);
                const node& blk = ensure_block(b);
                size_t i = block_start(b);
                return blk[p_idx - i];
            }

            virtual void str(size_t p_begin, size_t p_end, str_ranges& p_res) const
            {
                const size_t bb = block_num(p_begin);
                const size_t eb = block_num(p_end);
                for (size_t b = bb; b <= eb; ++b)
                {
                    size_t bs = block_start(b);
                    if (bs == p_end)
                    {
                        return;
                    }
                    const node& blk = ensure_block(b);
                    size_t bp = std::max(p_begin, bs);
                    size_t ep = std::min(p_end, block_start(b+1));
                    blk.str(bp - bs, ep - bs, p_res);
                }
            }

        private:
            const node& ensure_block(size_t p_block_num) const
            {
                if (m_blocks.find(p_block_num) == m_blocks.end())
                {
                    if (m_blocks.size() > 10)
                    {
                        m_blocks.clear();
                    }
                    m_blocks[p_block_num] = node_ptr(new atom(m_func(p_block_num)));
                }
                return *m_blocks[p_block_num];
            }

            size_t block_num(size_t p_pos) const
            {
                return p_pos >> m_block_bits;
            }

            size_t block_start(size_t p_block_num) const
            {
                return p_block_num << m_block_bits;
            }

            const size_t m_block_bits;
            std::function<std::string(size_t)> m_func;
            mutable std::unordered_map<size_t,node_ptr> m_blocks;
        };

        rope()
            : m_root(new atom(""))
        {
        }

        rope(const std::string& p_str)
            : m_root(new atom(p_str))
        {
        }

        rope(const rope& p_lhs, const rope& p_rhs)
            : m_root(new concat(p_lhs.m_root, p_rhs.m_root))
        {
        }

        rope(const rope& p_parent, size_t p_begin, size_t p_end)
            : m_root(new substr(p_parent.m_root, p_begin, p_end))
        {
        }

        rope(const size_t p_size, const size_t p_block_bits, std::function<std::string(size_t)> p_func)
            : m_root(new block(p_size, p_block_bits, p_func))
        {
        }

        size_t size() const
        {
            return m_root->size();
        }

        char operator[](size_t p_idx) const
        {
            return (*m_root)[p_idx];
        }

        std::string str() const
        {
            str_ranges rngs;
            size_t z = m_root->size();
            m_root->str(0, z, rngs);

            std::string res;
            res.reserve(z);
            for (auto i = rngs.begin(); i != rngs.end(); ++i)
            {
                res.insert(res.end(), i->first, i->second);
            }
            return res;
        }

        rope operator+(const rope& p_other) const
        {
            return rope(*this, p_other);
        }

        rope& operator +=(const rope& p_other)
        {
            m_root = node_ptr(new concat(m_root, p_other.m_root));
            return *this;
        }

        size_t find(char p_ch) const
        {
            size_t z = size();
            for (size_t i = 0; i < z; ++i)
            {
                if ((*m_root)[i] == p_ch)
                {
                    return i;
                }
            }
            return z;
        }

        rope slice(size_t p_begin, size_t p_end) const
        {
            return rope(*this, p_begin, p_end);
        }

        std::pair<rope,rope> split(size_t p_pos) const
        {
            size_t z = size();
            rope lhs = rope(*this, 0, p_pos);
            rope rhs = rope(*this, p_pos, z);
            return std::pair<rope,rope>(lhs, rhs);
        }

        void split(char p_ch, std::vector<rope>& p_parts) const
        {
            size_t z = size();
            size_t p = 0;
            for (size_t i = 0; i < z; ++i)
            {
                if ((*m_root)[i] == p_ch)
                {
                    p_parts.push_back(rope(*this, p, i));
                    p = i + 1;
                }
            }
            if (p != z)
            {
                p_parts.push_back(rope(*this, p, z));
            }
        }

    private:
        node_ptr m_root;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_ROPE_HPP
