#include "varoom/seq/genbank.hpp"

#include <initializer_list>
#include <fstream>
#include <sstream>
#include <iostream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE genbank tests
#include <boost/test/unit_test.hpp>

namespace // anonymous
{
    std::string compose(std::initializer_list<const char*> p_parts)
    {
        std::string s;
        for (auto i = p_parts.begin(); i != p_parts.end(); ++i)
        {
            s.insert(s.size(), *i);
        }
        return s;
    }
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( test_parser_1 )
{
    using namespace std;
    using namespace varoom::seq;

    const std::string fa = compose({
        "LOCUS       NC_000014          107349540 bp    DNA     linear   CON 13-AUG-2013\n",
        "DEFINITION  Homo sapiens chromosome 14, GRCh37.p13 Primary Assembly.\n",
        "ACCESSION   NC_000014\n",
        "VERSION     NC_000014.8\n",
        "//\n"
    });

    istringstream in(fa);
    size_t n = 0;
    for (genbank_reader r(in); r.more(); ++r, ++n)
    {
        switch (n)
        {
            case 0:
            {
                const genbank_record& rec = *r;
                BOOST_REQUIRE_EQUAL(rec.entries.size(), 4);
                BOOST_CHECK_EQUAL(rec.entries[0].name, "LOCUS");
                BOOST_CHECK_EQUAL(rec.entries[0].value, "NC_000014          107349540 bp    DNA     linear   CON 13-AUG-2013");
                BOOST_CHECK_EQUAL(rec.entries[1].name, "DEFINITION");
                BOOST_CHECK_EQUAL(rec.entries[1].value, "Homo sapiens chromosome 14, GRCh37.p13 Primary Assembly.");
                BOOST_CHECK_EQUAL(rec.entries[2].name, "ACCESSION");
                BOOST_CHECK_EQUAL(rec.entries[2].value, "NC_000014");
                BOOST_CHECK_EQUAL(rec.entries[3].name, "VERSION");
                BOOST_CHECK_EQUAL(rec.entries[3].value, "NC_000014.8");
                break;
            }
            default:
            {
                BOOST_REQUIRE_EQUAL(n == 0, true);
            }
        };
    }
}

BOOST_AUTO_TEST_CASE( test_parser_2 )
{
    using namespace std;
    using namespace varoom::seq;

    const std::string fa = compose({
        "SOURCE      Homo sapiens (human)\n",
        "  ORGANISM  Homo sapiens\n",
        "            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;\n",
        "            Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;\n",
        "            Catarrhini; Hominidae; Homo.\n",
        "//\n"
    });

    istringstream in(fa);
    size_t n = 0;
    for (genbank_reader r(in); r.more(); ++r, ++n)
    {
        switch (n)
        {
            case 0:
            {
                const genbank_record& rec = *r;
                BOOST_REQUIRE_EQUAL(rec.entries.size(), 1);
                BOOST_CHECK_EQUAL(rec.entries[0].name, "SOURCE");
                BOOST_CHECK_EQUAL(rec.entries[0].value, "Homo sapiens (human)");
                BOOST_REQUIRE_EQUAL(rec.entries[0].sub_entries.size(), 1);
                BOOST_REQUIRE_EQUAL(rec.entries[0].sub_entries[0].first, "ORGANISM");
                break;
            }
            default:
            {
                BOOST_REQUIRE_EQUAL(n == 0, true);
            }
        };
    }
}

BOOST_AUTO_TEST_CASE( test_parser_3 )
{
    using namespace std;
    using namespace varoom::seq;

    const std::string fa = compose({
        "REFERENCE   1  (bases 1 to 107349540)\n",
        "  CONSRTM   International Human Genome Sequencing Consortium\n",
        "  TITLE     Finishing the euchromatic sequence of the human genome\n",
        "  JOURNAL   Nature 431 (7011), 931-945 (2004)\n",
        "   PUBMED   15496913\n",
        "//\n"
    });

    istringstream in(fa);
    size_t n = 0;
    for (genbank_reader r(in); r.more(); ++r, ++n)
    {
        switch (n)
        {
            case 0:
            {
                const genbank_record& rec = *r;
                BOOST_REQUIRE_EQUAL(rec.entries.size(), 1);
                BOOST_CHECK_EQUAL(rec.entries[0].name, "REFERENCE");
                BOOST_CHECK_EQUAL(rec.entries[0].value, "1  (bases 1 to 107349540)");
                BOOST_REQUIRE_EQUAL(rec.entries[0].sub_entries.size(), 4);
                BOOST_REQUIRE_EQUAL(rec.entries[0].sub_entries[0].first, "CONSRTM");
                BOOST_REQUIRE_EQUAL(rec.entries[0].sub_entries[0].second, "International Human Genome Sequencing Consortium");
                BOOST_REQUIRE_EQUAL(rec.entries[0].sub_entries[3].first, "PUBMED");
                BOOST_REQUIRE_EQUAL(rec.entries[0].sub_entries[3].second, "15496913");
                break;
            }
            default:
            {
                BOOST_REQUIRE_EQUAL(n == 0, true);
            }
        };
    }
}

BOOST_AUTO_TEST_CASE( test_parser_4 )
{
    using namespace std;
    using namespace varoom::seq;

    const std::string fa = compose({
       //0123456789012345678901234567890123456789
        "FEATURES             Location/Qualifiers\n",
        "     source          107349540\n",
        "//\n"
    });

    istringstream in(fa);
    size_t n = 0;
    for (genbank_reader r(in); r.more(); ++r, ++n)
    {
        switch (n)
        {
            case 0:
            {
                const genbank_record& rec = *r;
                BOOST_REQUIRE_EQUAL(rec.entries.size(), 1);
                BOOST_CHECK_EQUAL(rec.entries[0].name, "FEATURES");
                BOOST_CHECK_EQUAL(rec.entries[0].features.size(), 1);
                BOOST_CHECK_EQUAL(rec.entries[0].features[0].type, "source");
                break;
            }
            default:
            {
                BOOST_REQUIRE_EQUAL(n == 0, true);
            }
        };
    }
}

BOOST_AUTO_TEST_CASE( test_parser_5 )
{
    using namespace std;
    using namespace varoom::seq;

    const std::string fa = compose({
       //0123456789012345678901234567890123456789
		"FEATURES             Location/Qualifiers\n",
		"     source          1..107349540\n",
		"                     /organism=\"Homo sapiens\"\n",
		"                     /mol_type=\"genomic DNA\"\n",
		"                     /db_xref=\"taxon:9606\"\n",
		"                     /chromosome=\"14\"\n",
        "//\n"
    });

    istringstream in(fa);
    size_t n = 0;
    for (genbank_reader r(in); r.more(); ++r, ++n)
    {
        switch (n)
        {
            case 0:
            {
                const genbank_record& rec = *r;
                BOOST_REQUIRE_EQUAL(rec.entries.size(), 1);
                BOOST_CHECK_EQUAL(rec.entries[0].name, "FEATURES");
                BOOST_CHECK_EQUAL(rec.entries[0].features.size(), 1);
                BOOST_CHECK_EQUAL(rec.entries[0].features[0].type, "source");
                BOOST_CHECK_EQUAL(rec.entries[0].features[0].qualifiers.size(), 4);
                BOOST_CHECK_EQUAL(rec.entries[0].features[0].qualifiers[0].first, "organism");
                BOOST_CHECK_EQUAL(rec.entries[0].features[0].qualifiers[0].second, "Homo sapiens");
                break;
            }
            default:
            {
                BOOST_REQUIRE_EQUAL(n == 0, true);
            }
        };
    }
}

BOOST_AUTO_TEST_CASE( test_parser_6 )
{
    using namespace std;
    using namespace varoom::seq;

    const std::string fa = compose({
       //0123456789012345678901234567890123456789
		"FEATURES             Location/Qualifiers\n",
		"     source          1..107349540\n",
		"     source          <1..107349540\n",
		"     source          1..>107349540\n",
		"     source          <1..>107349540\n",
        "//\n"
    });

    istringstream in(fa);
    size_t n = 0;
    for (genbank_reader r(in); r.more(); ++r, ++n)
    {
        switch (n)
        {
            case 0:
            {
                const genbank_record& rec = *r;
                BOOST_REQUIRE_EQUAL(rec.entries.size(), 1);
                BOOST_CHECK_EQUAL(rec.entries[0].name, "FEATURES");
                BOOST_CHECK_EQUAL(rec.entries[0].features.size(), 4);
                break;
            }
            default:
            {
                BOOST_REQUIRE_EQUAL(n == 0, true);
            }
        };
    }
}

BOOST_AUTO_TEST_CASE( test_parser_7 )
{
    using namespace std;
    using namespace varoom::seq;

    const std::string fa = compose({
       //0123456789012345678901234567890123456789
		"FEATURES             Location/Qualifiers\n",
		"     gene            complement(19410589..19413792)\n",
		"                     /gene=\"USP10P2\"\n",
		"                     /note=\"ubiquitin specific peptidase10 pseudogene 2;\n",
		"                     Derived by automated computational analysis using gene\n",
		"                     prediction method: Curated Genomic.\"\n",
		"                     /pseudo\n",
		"                     /db_xref=\"GeneID:780778\"\n",
		"                     /db_xref=\"HGNC:HGNC:33272\"\n",
        "//\n"
    });

    istringstream in(fa);
    size_t n = 0;
    for (genbank_reader r(in); r.more(); ++r, ++n)
    {
        switch (n)
        {
            case 0:
            {
                const genbank_record& rec = *r;
                BOOST_REQUIRE_EQUAL(rec.entries.size(), 1);
                BOOST_CHECK_EQUAL(rec.entries[0].name, "FEATURES");
                BOOST_CHECK_EQUAL(rec.entries[0].features.size(), 1);
                BOOST_CHECK_EQUAL(rec.entries[0].features[0].qualifiers.size(), 5);
                break;
            }
            default:
            {
                BOOST_REQUIRE_EQUAL(n == 0, true);
            }
        };
    }
}

BOOST_AUTO_TEST_CASE( test_parser_9 )
{
    using namespace std;
    using namespace varoom::seq;

    const std::string fa = compose({
		"ORIGIN      \n",
		"        1 cgttatttaa ggtgttacat agttctatgg aaatagggtc tatacctttc gccttacaat\n",
		"       61 gtaatttctt ..........\n",
        "//\n"
       //0123456789012345678901234567890123456789
    });

    istringstream in(fa);
    size_t n = 0;
    for (genbank_reader r(in); r.more(); ++r, ++n)
    {
        switch (n)
        {
            case 0:
            {
                const genbank_record& rec = *r;
                BOOST_REQUIRE_EQUAL(rec.entries.size(), 1);
                BOOST_CHECK_EQUAL(rec.entries[0].name, "ORIGIN");
                BOOST_CHECK_EQUAL(rec.entries[0].features.size(), 0);
                BOOST_CHECK_EQUAL(rec.entries[0].sequence, "cgttatttaaggtgttacatagttctatggaaatagggtctatacctttcgccttacaatgtaatttctt..........");
                break;
            }
            default:
            {
                BOOST_REQUIRE_EQUAL(n == 0, true);
            }
        };
    }
}

BOOST_AUTO_TEST_CASE( test_parser_10 )
{
    using namespace std;
    using namespace varoom::seq;

    ifstream in("tests/data/NC_001802.gbk");
    size_t n = 0;
    for (genbank_reader r(in); r.more(); ++r, ++n)
    {
        switch (n)
        {
            case 0:
            {
                const genbank_record& rec = *r;
                BOOST_REQUIRE_EQUAL(rec.entries.size(), 15);
                //cout << rec.to_json() << endl;
                break;
            }
            default:
            {
                BOOST_REQUIRE_EQUAL(n == 0, true);
            }
        };
    }
}
