#ifndef VAROOM_UTIL_LRU_CACHE_HPP
#define VAROOM_UTIL_LRU_CACHE_HPP

namespace varoom
{
    template <typename K, typename V>
    class lru_cache
    {
    public:
        lru_cache(const size_t& p_max_items,
                  std::function<V(K const&)> p_lookup)
            : m_max_items(p_max_items), m_lookup(p_lookup), m_seq_num(0)
        {
        }

        const V& operator[](const K& p_key)
        {
            auto itr = m_cache.find(p_key);
            if (itr == m_cache.end())
            {
                load(p_key);
                itr = m_cache.find(p_key);
            }
            touch(p_key);
            return itr->second;
        }

    private:

        void load(const K& p_key)
        {
            if (m_cache.size() == m_max_items)
            {
                evict();
            }
            m_cache[p_key] = m_lookup(p_key);
        }

        void touch(p_key)
        {
            auto itr = m_lru_idx.find(p_key);
            if (itr != m_lru_idx.end())
            {
                m_lru.erase(itr->second);
            }
            itr->second = ++m_seq_num;
            m_lru[itr->second] = itr->first;
        }

        void evict()
        {
            auto itr = m_lru.begin();
            if (itr == m_lru.end())
            {
                return;
            }
            m_lru_idx.erase(itr->second);
            m_cache.erase(itr->second);
            m_lru.erase(itr);
        }

        const size_t m_max_items;
        std::function<V(K const&)> m_lookup;
        size_t m_seq_num;
        std::map<size_t,K> m_lru;
        std::map<K,size_t> m_lru_idx;
        std::map<K,V> m_cache;
    };
}
// namespace varoom

#endif // VAROOM_UTIL_LRU_CACHE_HPP
