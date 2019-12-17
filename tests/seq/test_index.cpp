#include "varoom/seq/index.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE seq_index tests
#include <boost/test/unit_test.hpp>

#include <random>

namespace // anonymous
{
    std::map<std::string,std::string>& fs()
    {
        static std::map<std::string,std::string> m;
        if (m.size() == 0)
        {
            
            std::mt19937_64 rng(17);

            std::uniform_int_distribution<size_t> Z(64*1024, 128*1024);
            std::uniform_int_distribution<int> U(0, 4);

            std::vector<std::string> chrs = {"chr1", "chr2", "chr3", "chr4"};

            std::ostringstream toc;

            for (size_t i = 0; i < chrs.size(); ++i)
            {
                size_t z = Z(rng);

                std::cerr << chrs[i] << '\t' << z << std::endl;

                std::string x;
                x.reserve(z);

                for (size_t j = 0; j < z; ++j)
                {
                    x.push_back("ACGTN"[U(rng)]);
                }
                std::string fn = std::string("/varoom/hg19/") + chrs[i] + std::string(".txt");
                m[fn] = x;
                toc << chrs[i] << '\t' << z << '\t' << fn << std::endl;
            }
            m["/varoom/hg19/toc.txt"] = toc.str();
            std::cerr << m["/varoom/hg19/toc.txt"];
        }
        return m;
    }
}
// namespace anonymous

BOOST_AUTO_TEST_CASE( enc1 )
{
    const std::string& s = fs().find("/varoom/hg19/chr2.txt")->second;
    varoom::files::string_impl(fs());
    varoom::seq::index I("hg19/toc.txt", 10);

    {
        std::string x(s.begin() + 1234, s.begin() + 5678);
        std::string r = I.get("chr2", 1234, 5678);
        BOOST_CHECK_EQUAL(r, x);
    }
}
