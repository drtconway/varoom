#ifndef VAROOM_SEQ_COMPACT_SEQ_HPP
#define VAROOM_SEQ_COMPACT_SEQ_HPP

#ifndef VAROOM_UTIL_FILES_HPP
#include "varoom/util/files.hpp"
#endif

#ifndef VAROOM_UTIL_GZIP_STRINGS_HPP
#include "varoom/util/gzip_strings.hpp"
#endif

#ifndef VAROOM_SEQ_FASTA_HPP
#include "varoom/seq/fasta.hpp"
#endif

#include <nlohmann/json.hpp>

namespace varoom
{
    namespace seq
    {
        class compact_seq
        {
        public:
            typedef nlohmann::json json;
            typedef std::map<std::string,std::vector<std::pair<size_t,size_t>>> toc_type;

            static const size_t B = 1024*1024;

            compact_seq(const std::string& p_src_base_name)
                : m_dat_holder(files::in(p_src_base_name + ".dat")), m_dat(**m_dat_holder)
            {
                input_file_holder_ptr tocp = files::in(p_src_base_name + ".toc.gz");
                json toc;
                (**tocp) >> toc;
                for (json::iterator itr = toc.begin(); itr != toc.end(); ++itr)
                {
                    std::string acc = itr.key();
                    const json& ranges = itr.value();
                    for (size_t i = 0; i < ranges.size(); ++i)
                    {
                        std::pair<size_t,size_t> v(ranges[i][0], ranges[i][1]);
                        m_toc[acc].push_back(v);
                    }
                }
            }

            static void make(const std::string& p_src_name, const std::string& p_dest_base_name)
            {
                toc_type toc;
                input_file_holder_ptr inp = files::in(p_src_name);
                output_file_holder_ptr outp = files::out(p_dest_base_name + ".dat");
                for (fasta_reader r(**inp); r.more(); ++r)
                {
                    const fasta_read& rd = *r;
                    size_t ztot = 0;
                    for (size_t p = 0; p < rd.second.size(); p += B)
                    {
                        size_t q = std::min(rd.second.size(), p + B);
                        std::string s(rd.second.begin() + p, rd.second.begin() + q);
                        std::string z = gzip_strings::compress(s);
                        (**outp).write(&z[0], z.size());
                        ztot += z.size();
                        toc[rd.first].push_back(std::pair<size_t,size_t>(ztot, ztot+z.size()));
                    }
                }

                json res = json::object();
                for (auto itr = toc.begin(); itr != toc.end(); ++itr)
                {
                    json vec = json::array();
                    for (auto jtr = itr->second.begin(); jtr != itr->second.end(); ++jtr)
                    {
                        vec.push_back(json{jtr->first, jtr->second});
                    }
                    res[itr->first] = vec;
                }
                {
                    output_file_holder_ptr tocp = files::out(p_dest_base_name + ".toc.gz");
                    (**tocp) << res;
                }
            }

            std::string get(const std::string& p_acc, const size_t& p_begin, const size_t& p_end)
            {
                size_t b1 = block_of(p_begin);
                size_t b2 = block_of(p_end);
                if (b1 == b2)
                {
                    std::string blk;
                    load_block(p_acc, b1, blk);
                    std::pair<size_t,size_t> r = block_range(b1);
                    size_t b = p_begin - r.first;
                    size_t e = p_end - r.first;
                    return std::string(blk.begin() + b, blk.begin() + e);
                }
                
                std::string blk;
                load_block(p_acc, b1, blk);
                std::pair<size_t,size_t> r1 = block_range(b1);
                size_t b = p_begin - r1.first;
                std::string res(blk.begin() + b, blk.end());
                for (size_t bn = b1 + 1; bn < b2; ++bn)
                {
                    load_block(p_acc, bn, blk);
                    res.insert(res.end(), blk.begin(), blk.end());
                }
                load_block(p_acc, b2, blk);
                std::pair<size_t,size_t> r2 = block_range(b2);
                size_t e = p_end - r2.first;
                res.insert(res.end(), blk.begin(), blk.begin() + e);
                return res;
            }

        private:
            static size_t block_of(const size_t& p_position)
            {
                return p_position / B;
            }

            static std::pair<size_t,size_t> block_range(const size_t& p_block_num)
            {
                size_t b = B * p_block_num;
                size_t e = B * (p_block_num + 1);
                return std::pair<size_t,size_t>(b, e);
            }

            void load_block(const std::string& p_acc, const size_t& p_block_num, std::string& p_res)
            {
                std::pair<size_t,size_t> rng = m_toc.find(p_acc)->second[p_block_num];
                size_t l = rng.second - rng.first;
                std::string buf(l, '\0');
                m_dat.seekg(rng.first);
                m_dat.read(&buf[0], l);
                p_res = gzip_strings::decompress(buf);
            }

            input_file_holder_ptr m_dat_holder;
            std::istream& m_dat;
            toc_type m_toc;
        };
    }
    // namespace seq
}
// namespace varoom

#endif // VAROOM_SEQ_COMPACT_SEQ_HPP
