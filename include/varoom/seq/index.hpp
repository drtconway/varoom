#ifndef VAROOM_SEQ_INDEX_HPP
#define VAROOM_SEQ_INDEX_HPP

#include <string>
#include <vector>

#ifndef VAROOM_UTIL_FILES_HPP
#include "varoom/util/files.hpp"
#endif

#ifndef VAROOM_UTIL_ROPE_HPP
#include "varoom/util/rope.hpp"
#endif

#ifndef VAROOM_UTIL_TABLE_HPP
#include "varoom/util/table.hpp"
#endif

namespace varoom
{
    namespace seq
    {
        namespace detail
        {
            struct genome : varoom::table<std::string,uint64_t,std::string>
            {
                using tuple = tuple_type;
                using table = basic_inmemory_table;
                using istream_reader = basic_istream_table;
                using ostream_writer = basic_ostream_table;
                using read_iterator = basic_read_iterator;
                using write_iterator = basic_write_iterator;

                static constexpr std::size_t chr = 0;
                static constexpr std::size_t len = 1;
                static constexpr std::size_t pth = 2;
            };
        }
        // namespace detail

        class index
        {
        public:
            index(const std::string& p_toc_filename, size_t p_block_size_bits)
                : m_block_size_bits(p_block_size_bits)
            {
                boost::filesystem::path toc_path = boost::filesystem::absolute(p_toc_filename);
                m_parent_path = toc_path.parent_path();

                std::function<void(const detail::genome::tuple&)> f = [this] (const detail::genome::tuple& p_item) mutable {
                    std::string chr = std::get<detail::genome::chr>(p_item);
                    uint64_t len = std::get<detail::genome::len>(p_item);
                    std::string pth = std::get<detail::genome::pth>(p_item);
                    boost::filesystem::path p(pth);
                    boost::filesystem::path r = p;
                    if (r.is_relative())
                    {
                        r = m_parent_path;
                        r /= p;
                    }

                    size_t B = m_block_size_bits;
                    std::string fn = r.string();

                    std::function<std::string(size_t)> ldr = [B,fn](size_t n) {
                        input_file_holder_ptr inp = files::in(fn);
                        std::istream& in = **inp;

                        std::string res;
                        size_t z = 1ULL << B;
                        res.resize(z);

                        size_t pos = n << B;
                        in.seekg(pos);
                        in.read(&res[0], z);
                        return res;
                    };

                    m_index[chr] = varoom::rope(len, B, ldr);
                };
                input_file_holder_ptr tocp = files::in(toc_path.string());
                detail::genome::istream_reader g(**tocp, true);
                table_utils::for_each(g, f);
            }

            std::string get(const std::string& p_acc, size_t p_begin, size_t p_end)
            {
                auto itr = m_index.find(p_acc);
                if (itr == m_index.end())
                {
                    throw std::runtime_error("accession not found");
                }
                return itr->second.slice(p_begin, p_end).str();
            }

        private:
            const size_t m_block_size_bits;
            boost::filesystem::path m_parent_path;
            std::map<std::string,varoom::rope> m_index;
        };
    }
    // namespace seq
}
// namespace varoom

#endif // VAROOM_SEQ_INDEX_HPP
