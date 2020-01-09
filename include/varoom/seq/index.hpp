#ifndef VAROOM_SEQ_INDEX_HPP
#define VAROOM_SEQ_INDEX_HPP

#include <string>
#include <vector>

#ifndef VAROOM_UTIL_FILES_HPP
#include "varoom/util/files.hpp"
#endif

#ifndef VAROOM_SEQ_FASTA_HPP
#include "varoom/seq/fasta.hpp"
#endif

#ifndef VAROOM_SEQ_KIMMO_HPP
#include "varoom/seq/kimmo.hpp"
#endif

namespace varoom
{
    namespace seq
    {
        namespace detail
        {
        }
        // namespace detail

        class index
        {
        public:
            using toc_type = std::unordered_map<std::string,size_t>;

            index(const std::string& p_name)
                : m_inp(varoom::files::in(p_name)), m_in(**m_inp)
            {
                load();
            }

            index(std::istream& p_in)
                : m_in(p_in)
            {
                load();
            }

            std::string get(const std::string& p_acc, size_t p_begin, size_t p_end)
            {
                load_accession(p_acc);
                return m_curr_seq->slice(p_begin, p_end);
            }

            void load()
            {
                m_in.seekg(0, std::ios_base::end);
                uint64_t l = m_in.tellg();
                m_in.seekg(l - sizeof(uint64_t), std::ios_base::beg);
                uint64_t z;
                m_in.read(reinterpret_cast<char*>(&z), sizeof(uint64_t));
                m_in.seekg(z);
                blob::load(m_in, [this](std::istream& i) mutable {
                    nlohmann::json jj;
                    i >> jj;
                    for (nlohmann::json::iterator itr = jj.begin(); itr != jj.end(); ++itr)
                    {
                        size_t x = itr.value();
                        m_toc[itr.key()] = x;
                    }
                });
            }

            void load_accession(const std::string& p_acc)
            {
                if (m_curr_acc == p_acc && m_curr_seq)
                {
                    return;
                }
                auto itr = m_toc.find(p_acc);
                if (itr == m_toc.end())
                {
                    throw std::runtime_error("accession not found: " + p_acc);
                }
                uint64_t pos = itr->second;
                m_in.seekg(pos, std::ios_base::beg);
                m_curr_acc = p_acc;
                m_curr_seq = std::shared_ptr<varoom::seq::kimmo>(new kimmo(m_in));
            }

            static void make(const std::vector<std::string>& p_sources, const std::string& p_outname)
            {
                varoom::output_file_holder_ptr outp = varoom::files::out(p_outname);
                make(p_sources, **outp);
            }

            static void make(const std::vector<std::string>& p_sources, std::ostream& p_out)
            {
                toc_type toc;
                for (size_t i = 0; i < p_sources.size(); ++i)
                {
                    std::cerr << "processing " << p_sources[i] << " @ " << p_out.tellp() << std::endl;
                    input_file_holder_ptr inp = files::in(p_sources[i]);
                    for (varoom::seq::fasta_reader r(**inp); r.more(); ++r)
                    {
                        const fasta_read rd = *r;
                        size_t off = p_out.tellp();
                        toc[rd.first] = off;
                        kimmo::make(rd.second, p_out);
                    }
                }
                uint64_t off = p_out.tellp();
                blob::save(p_out, [&](std::ostream& o) {
                    nlohmann::json jj = toc;
                    o << jj;
                });
                p_out.write(reinterpret_cast<const char*>(&off), sizeof(off));
            }

        private:
            varoom::input_file_holder_ptr m_inp;
            std::istream& m_in;
            toc_type m_toc;
            std::string m_curr_acc;
            std::shared_ptr<varoom::seq::kimmo> m_curr_seq;
        };
    }
    // namespace seq
}
// namespace varoom

#endif // VAROOM_SEQ_INDEX_HPP
