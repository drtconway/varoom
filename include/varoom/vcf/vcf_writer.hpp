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

// ##FILTER=<ID=PASS,Description="All filters passed">
// ##FILTER=<ID=LowQual,Description="Low quality">
// ##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
// ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
// ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
// ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
// ##FORMAT=<ID=PMCAD,Number=A,Type=Integer,Description="Depth of alternate-supporting bases - Calculated By Bioinformatics Dept">
// ##FORMAT=<ID=PMCADF,Number=A,Type=String,Description="Depth of alternate-supporting bases on forward strand - Calculated By Bioinformatics Dept">
// ##FORMAT=<ID=PMCADR,Number=A,Type=String,Description="Depth of alternate-supporting bases on reverse strand - Calculated By Bioinformatics Dept">
// ##FORMAT=<ID=PMCBDIR,Number=A,Type=String,Description="T/F Indicating if variant is bidirectional (N/A if no alt reads) - Calculated By Bioinformatics Dept">
// ##FORMAT=<ID=PMCDP,Number=1,Type=Integer,Description="Total read depth (includes bases supporting other alleles) - Calculated By Bioinformatics Dept">
// ##FORMAT=<ID=PMCFREQ,Number=A,Type=Float,Description="Variant allele frequency - Calculated By Bioinformatics Dept">
// ##FORMAT=<ID=PMCRD,Number=1,Type=Integer,Description="Depth of reference-supporting bases - Calculated By Bioinformatics Dept">
// ##FORMAT=<ID=PMCRDF,Number=1,Type=String,Description="Depth of reference-supporting bases on forward strand - Calculated By Bioinformatics Dept">
// ##FORMAT=<ID=PMCRDR,Number=1,Type=String,Description="Depth of reference-supporting bases on reverse strand - Calculated By Bioinformatics Dept">
// ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
// ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
// ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
// ##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
// ##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
// ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
// ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
// ##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
// ##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
// ##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
// ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
// ##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
// ##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
// ##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
// ##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
// ##reference=file:///data/seqliner/seqliner-resources/reference_genomes/g1k_v37/human_g1k_v37.fasta
// ##INFO=<ID=OLD_MULTIALLELIC,Number=1,Type=String,Description="Original chr:pos:ref:alt encoding">
// ##INFO=<ID=OLD_VARIANT,Number=.,Type=String,Description="Original chr:pos:ref:alt encoding">
// ##INFO=<ID=HGVSg,Number=1,Type=String,Description="HGVSg format of variant">
// ##INFO=<ID=HGVSc,Number=1,Type=String,Description="HGVSc format of variant">
// ##INFO=<ID=HGVSp,Number=1,Type=String,Description="HGVSp format of variant">
// ##INFO=<ID=gene,Number=1,Type=String,Description="gene of variant">
// ##INFO=<ID=lrg,Number=1,Type=String,Description="LRG transcript of variant">
// ##INFO=<ID=muterr,Number=1,Type=String,Description="Mutalyzer error message for variant">
// ##INFO=<ID=status,Number=1,Type=String,Description="Mutalyzer status of variant">
// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  10458703

namespace varoom
{
    namespace vcf
    {
        class vcf_writer : public vcf_handler
        {
        public:
            vcf_writer(std::ostream& p_out, const std::string& p_version = "VCF4.2")
                : m_out(p_out)
            {
                m_out << "##fileformat=" << p_version << std::endl;

            }

            vcf_writer(std::ostream& p_out, const std::string& p_version, const nlohmann::json& p_meta)
                : m_out(p_out)
            {
                m_out << "##fileformat=" << p_version << std::endl;
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
            }

            virtual void operator()(const std::string& p_chr,
                                    const std::int64_t& p_pos,
                                    const std::string& p_id,
                                    const std::string& p_ref,
                                    const std::string& p_alt,
                                    const std::int64_t& p_qual,
                                    const std::string& p_filter,
                                    const lazy_vcf_info& p_info,
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
                    for (auto itr = ifo.begin(); itr != ifo.end(); ++itr)
                    {
                        m_out << (itr == ifo.begin() ? '\t' : ';');
                        m_out << itr->first << '=' << itr->second;
                    }
                }
                const std::vector<vcf_info>& genotypes = p_genotypes.get();
                if (genotypes.size() > 0)
                {
                    const vcf_info& ifo = genotypes[0];
                    for (auto itr = ifo.begin(); itr != ifo.end(); ++itr)
                    {
                        m_out << (itr == ifo.begin() ? '\t' : ':');
                        m_out << itr->first;
                    }
                }
                for (size_t i = 0; i < genotypes.size(); ++i)
                {
                    const vcf_info& ifo = genotypes[i];
                    for (auto itr = ifo.begin(); itr != ifo.end(); ++itr)
                    {
                        m_out << (itr == ifo.begin() ? '\t' : ':');
                        m_out << itr->second;
                    }
                }
                m_out << std::endl;
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
        };
    }
    // namespace vcf
}
// namespace varoom


#endif // VAROOM_VCF_VCF_WRITER_HPP
