#ifndef VAROOM_UTIL_RANK_SET_HPP
#define VAROOM_UTIL_RANK_SET_HPP

#include <cassert>
#include <functional>
#include <memory>

#include <iostream>

namespace varoom
{
    namespace detail
    {
        template <typename KeyType>
        struct avl_tree;
        template <typename KeyType>
        using avl_tree_ptr = std::shared_ptr<avl_tree<KeyType>>;

        template <typename KeyType>
        struct avl_tree
        {
            using ptr_type = avl_tree_ptr<KeyType>;

            KeyType m_key;
            size_t m_size;
            int m_height;
            ptr_type m_lhs;
            ptr_type m_rhs;

            avl_tree(const KeyType& p_key)
                : m_key(p_key), m_size(1), m_height(1)
            {
            }

            avl_tree(const KeyType& p_key, ptr_type p_lhs, ptr_type p_rhs)
                : m_key(p_key), m_lhs(p_lhs), m_rhs(p_rhs)
            {
                refresh_size_and_height();
            }

            size_t size() const
            {
                return m_size;
            }

            size_t height() const
            {
                return m_height;
            }

            bool contains(const KeyType& p_key) const
            {
                if (m_key == p_key)
                {
                    return true;
                }
                else if (p_key < m_key)
                {
                    if (m_lhs)
                    {
                        return m_lhs->contains(p_key);
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    if (m_rhs)
                    {
                        return m_rhs->contains(p_key);
                    }
                    else
                    {
                        return false;
                    }
                }
            }

            size_t rank(const KeyType& p_key) const
            {
                if (p_key == m_key)
                {
                    return (m_lhs ? m_lhs->size() : 0);
                }
                if (p_key < m_key)
                {
                    return (m_lhs ? m_lhs->rank(p_key) : 0);
                }
                else
                {
                    return (m_lhs ? m_lhs->size() : 0) + 1 + (m_rhs ? m_rhs->rank(p_key) : 0);
                }
            }

            const KeyType& select(size_t p_pos) const
            {
                size_t lz = (m_lhs ? m_lhs->size() : 0);
                if (p_pos < lz)
                {
                    return m_lhs->select(p_pos);
                }
                p_pos -= lz;
                if (p_pos == 0)
                {
                    return m_key;
                }
                else
                {
                    return m_rhs->select(p_pos - 1);
                }
            }

            void refresh_size_and_height()
            {
                m_size = 1;
                int lh = 0;
                int rh = 0;
                if (m_lhs)
                {
                    m_size += m_lhs->size();
                    lh = m_lhs->height();
                }
                if (m_rhs)
                {
                    m_size += m_rhs->size();
                    rh = m_rhs->height();
                }
                m_height = 1 + std::max(lh, rh);
            }

            int balance_factor() const
            {
                int lh = 0;
                int rh = 0;
                if (m_lhs)
                {
                    lh = m_lhs->height();
                }
                if (m_rhs)
                {
                    rh = m_rhs->height();
                }
                return lh - rh;
            }

            static ptr_type rotate_left(ptr_type x)
            {
                assert(x);
                assert(x->m_rhs);
                ptr_type z = x->m_rhs;
                ptr_type zl = z->m_lhs;
                x->m_rhs = zl;
                x->refresh_size_and_height();
                z->m_lhs = x;
                z->refresh_size_and_height();
                return z;
            }

            static ptr_type rotate_right(ptr_type x)
            {
                assert(x);
                assert(x->m_lhs);
                ptr_type z = x->m_lhs;
                ptr_type zr = z->m_rhs;
                x->m_lhs = zr;
                x->refresh_size_and_height();
                z->m_rhs = x;
                z->refresh_size_and_height();
                return z;
            }

            static ptr_type insert(const KeyType& p_key, ptr_type x)
            {
                if (!x)
                {
                    return ptr_type(new avl_tree(p_key));
                }

                if (x->m_key == p_key)
                {
                    return x;
                }

                if (p_key < x->m_key)
                {
                    x->m_lhs = insert(p_key, x->m_lhs);
                }
                else
                {
                    x->m_rhs = insert(p_key, x->m_rhs);
                }

                x->refresh_size_and_height();

                int b = x->balance_factor();

                assert(-2 <= b && b <= 2);

                if (b > 1 && p_key < x->m_lhs->m_key)
                {
                    // left-left
                    return rotate_right(x);
                }
                if (b < -1 && x->m_rhs->m_key < p_key)
                {
                    // right-right
                    return rotate_left(x);
                }
                if (b > 1 && x->m_lhs->m_key < p_key)
                {
                    // left-right
                    x->m_lhs = rotate_left(x->m_lhs);
                    return rotate_right(x);
                }
                if (b < -1 && p_key < x->m_rhs->m_key)
                {
                    // right-left
                    x->m_rhs = rotate_right(x->m_rhs);
                    return rotate_left(x);
                }
                return x;
            }

            static ptr_type min_value_node(ptr_type x)
            {
                while (x->m_lhs)
                {
                    x = x->m_lhs;
                }
                return x;
            }

            static ptr_type remove(const KeyType& p_key, ptr_type x)
            {
                if (!x)
                {
                    // Key wasn't there.
                    return x;
                }

                if (p_key == x->m_key)
                {
                    if (!x->m_lhs || !x->m_rhs)
                    {
                        return (x->m_lhs ? x->m_lhs : x->m_rhs);
                    }
                    ptr_type t = min_value_node(x->m_rhs);
                    x->m_key = t->m_key;
                    x->m_rhs = remove(x->m_key, x->m_rhs);
                }
                else if (p_key < x->m_key)
                {
                    x->m_lhs = remove(p_key, x->m_lhs);
                }
                else
                {
                    x->m_rhs = remove(p_key, x->m_rhs);
                }

                x->refresh_size_and_height();
                int b = x->balance_factor();

                if (b > 1)
                {
                    int lb = (x->m_lhs ? x->m_lhs->balance_factor() : 0);
                    if (lb < 0)
                    {
                        x->m_lhs = rotate_left(x->m_lhs);
                    }
                    return rotate_right(x);
                }
                if (b < -1)
                {
                    int rb = (x->m_rhs ? x->m_rhs->balance_factor() : 0);
                    if (rb > 0)
                    {
                        x->m_rhs = rotate_right(x->m_rhs);
                    }
                    return rotate_left(x);
                }
                return x;
            }

            static void sane(ptr_type x)
            {
                if (!x)
                {
                    return;
                }
                sane(x->m_lhs);
                sane(x->m_rhs);
                //int b = x->balance_factor();
                //assert(-1 <= b && b <= 1);
            }

            static void visit(ptr_type x, std::function<void(const KeyType&)> p_vis)
            {
                if (!x)
                {
                    return;
                }
                visit(x->m_lhs, p_vis);
                p_vis(x->m_key);
                visit(x->m_rhs, p_vis);
            }
        };
    }
    // namespace detail

    template <typename T>
    class rank_set
    {
    public:
        rank_set()
        {
        }

        size_t size() const
        {
            return (m_root ? m_root->size() : 0);
        }

        size_t rank(const T& p_item) const
        {
            return (m_root ? m_root->rank(p_item) : 0);
        }

        const T& select(size_t p_index) const
        {
            assert(p_index < size());
            return m_root->select(p_index);
        }

        rank_set& insert(const T& p_item)
        {
            m_root = varoom::detail::avl_tree<T>::insert(p_item, m_root);
            return *this;
        }

        rank_set& erase(const T& p_item)
        {
            m_root = varoom::detail::avl_tree<T>::remove(p_item, m_root);
            return *this;
        }

    private:
        varoom::detail::avl_tree_ptr<T> m_root;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_RANK_SET_HPP
