#ifndef VAROOM_SEQ_SEQUENCE_FACTORY_HPP
#define VAROOM_SEQ_SEQUENCE_FACTORY_HPP

#ifndef VAROOM_SEQ_FASTA_HPP
#include "varoom/seq/fasta.hpp"
#endif

#ifndef VAROOM_UTIL_FILES_HPP
#include "varoom/util/files.hpp"
#endif

namespace varoom
{
    namespace seq
    {
        class sequence_factory
        {
        public:
            sequence_factory(const std::string& p_genome_directory)
                : m_genome_directory(p_genome_directory)
            {
            }

            const fasta_read& operator[](const std::string& p_accession)
            {
                if (m_curr.first != p_accession)
                {
                    std::string path = get_path(p_accession);
                    load_accession(path);
                }
                return m_curr;
            }

        private:
            std::string get_path(const std::string& p_accession) const
            {
                std::string path = m_genome_directory + "/" + p_accession + ".fa.gz";
                if (files::exists(path))
                {
                    return path;
                }
                path = m_genome_directory + "/chr" + p_accession + ".fa.gz";
                if (files::exists(path))
                {
                    return path;
                }
                path = m_genome_directory + "/" + p_accession + ".fa";
                if (files::exists(path))
                {
                    return path;
                }
                path = m_genome_directory + "/chr" + p_accession + ".fa";
                if (files::exists(path))
                {
                    return path;
                }
                throw std::runtime_error("could not locate sequence data");
            }

            void load_accession(const std::string& p_accession_name)
            {
                input_file_holder_ptr inp = files::in(p_accession_name);
                for (fasta_reader r(**inp); r.more(); ++r)
                {
                    m_curr = *r;
                    return;
                }
                throw std::runtime_error("no sequence found");
            }

            const std::string m_genome_directory;
            fasta_read m_curr;
        };
    }
    // namespace seq
}
// namespace varoom

#endif // VAROOM_SEQ_SEQUENCE_FACTORY_HPP
