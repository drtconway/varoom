#ifndef VAROOM_VCF_VCF_WRITER_HPP
#define VAROOM_VCF_VCF_WRITER_HPP

#ifndef VAROOM_VCF_VCF_HANDLER_HPP
#include "varoom/vcf/vcf_handler.hpp"
#endif

#ifndef VAROOM_UTIL_TEXT_HPP
#include "varoom/util/text.hpp"
#endif

#include <initializer_list>
#include <ostream>
#include <nlohmann/json.hpp>

namespace varoom
{
    namespace vcf
    {
        class vcf_writer : public vcf_handler
        {
        public:
            vcf_writer(std::ostream& p_out, const std::string& p_version = "VCF4.2")
                : m_out(p_out), m_no_more_meta(false)
            {
                m_out << "##fileformat=" << p_version << std::endl;

            }

            void meta(const nlohmann::json& p_meta)
            {
                if (m_no_more_meta)
                {
                    throw std::runtime_error("no more meta data can be written.");
                }
                if (p_meta.count("FILTER"))
                {
                    write_filters(p_meta["FILTER"]);
                }
                if (p_meta.count("INFO"))
                {
                    write_infos(p_meta["INFO"]);
                }
                if (p_meta.count("FORMAT"))
                {
                    write_formats(p_meta["FORMAT"]);
                }
            }

            void samples(const std::vector<std::string>& p_sample_names)
            {
                m_out << "#" << varoom::text::tabs({"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"});
                if (p_sample_names.size() > 0)
                {
                    m_out << '\t' << "FORMAT";
                }
                for (size_t i = 0; i < p_sample_names.size(); ++i)
                {
                    m_out << '\t' << p_sample_names[i];
                }
                m_out << std::endl;
                m_no_more_meta = true;
            }

            virtual void operator()(const std::string& p_chr,
                                    const std::int64_t& p_pos,
                                    const std::string& p_id,
                                    const std::string& p_ref,
                                    const std::string& p_alt,
                                    const double& p_qual,
                                    const std::string& p_filter,
                                    const lazy<vcf_info>& p_info,
                                    const lazy<std::vector<vcf_info>>& p_genotypes)
            {
                m_out << p_chr
                      << '\t' << p_pos
                      << '\t' << p_id
                      << '\t' << p_ref
                      << '\t' << p_alt
                      << '\t' << p_qual
                      << '\t' << p_filter;
                {
                    const vcf_info& ifo = p_info.get();
                    for (size_t i = 0; i < ifo.size(); ++i)
                    {
                        m_out << (i == 0 ? '\t' : ';');
                        m_out << ifo[i].first << '=' << ifo[i].second;
                    }
                }
                const std::vector<vcf_info>& genotypes = p_genotypes.get();
                if (genotypes.size() > 0)
                {
                    const vcf_info& ifo = genotypes[0];
                    for (size_t i = 0; i < ifo.size(); ++i)
                    {
                        m_out << (i == 0 ? '\t' : ':');
                        m_out << ifo[i].first;
                    }
                }
                for (size_t i = 0; i < genotypes.size(); ++i)
                {
                    const vcf_info& ifo = genotypes[i];
                    for (size_t i = 0; i < ifo.size(); ++i)
                    {
                        m_out << (i == 0 ? '\t' : ':');
                        m_out << ifo[i].second;
                    }
                }
                m_out << std::endl;
                m_no_more_meta = true;
            }

            virtual void error(const size_t& p_line_no, const std::string& p_line, const std::string& p_message)
            {
            }

        private:
            void write_filters(const nlohmann::json& p_filters)
            {
                for (size_t i = 0; i < p_filters.size(); ++i)
                {
                    const nlohmann::json& itm = p_filters[i];
                    m_out << "##FILTER=<ID=" << itm["ID"] << ",Description=\"" << itm["Description"] << "\">"
                          << std::endl;
                }
            }

            void write_infos(const nlohmann::json& p_infos)
            {
                for (size_t i = 0; i < p_infos.size(); ++i)
                {
                    const nlohmann::json& itm = p_infos[i];
                    m_out << "##INFO=<"
                                << "ID=" << itm["ID"]
                                << ",Number=" << itm["Number"]
                                << ",Type=" << itm["Type"]
                                << ",Description=\"" << itm["Description"] << "\"";
                    if (itm.count("Source"))
                    {
                        m_out << ",Source=" << itm["Source"];
                    }
                    if (itm.count("Version"))
                    {
                        m_out << ",Version=" << itm["Version"];
                    }
                    m_out << ">" << std::endl;
                }
            }

            void write_formats(const nlohmann::json& p_formats)
            {
                for (size_t i = 0; i < p_formats.size(); ++i)
                {
                    const nlohmann::json& itm = p_formats[i];
                    m_out << "##FORMAT=<"
                                << "ID=" << itm["ID"]
                                << ",Number=" << itm["Number"]
                                << ",Type=" << itm["Type"]
                                << ",Description=\"" << itm["Description"] << "\">";
                    m_out << ">" << std::endl;
                }
            }

            std::ostream& m_out;
            bool m_no_more_meta;
        };
    }
    // namespace vcf
}
// namespace varoom


#endif // VAROOM_VCF_VCF_WRITER_HPP
