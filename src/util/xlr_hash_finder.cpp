#include "varoom/util/xlr_hash_finder.hpp"

#include <fstream>
#include <boost/algorithm/string/trim.hpp>

#define BOOST_LOG_DYN_LINK 1
#include <boost/log/trivial.hpp>

using namespace std;
using namespace boost;
using namespace varoom;

namespace // anonymous
{
    void read_strings(const string& p_file_name, vector<string>& p_strings)
    {
        BOOST_LOG_TRIVIAL(info) << "reading strings from " << p_file_name;
        ifstream in(p_file_name);
        string l;
        while (getline(in, l))
        {
            boost::algorithm::trim(l);
            p_strings.push_back(l);
        }
        BOOST_LOG_TRIVIAL(info) << p_strings.size() << " strings read.";
    }
}
// namespace anonymous

int main(int argc, const char** argv)
{
    vector<string> s;
    read_strings(argv[1], s);
    xlr_hash_finder f(2, 19, s);
    vector<uint64_t> coeffs;
    f.find(1000, coeffs);
    cout << varoom::detail::as_json(coeffs) << endl;
    return 0;
}
