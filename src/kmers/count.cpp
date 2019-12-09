#include "varoom/kmers.hpp"
#include "varoom/seq/fastq.hpp"
#include "varoom/seq/fastq_pair.hpp"
#include "varoom/util/files.hpp"
#include "varoom/command.hpp"

#include <chrono>
#include <unordered_map>

int main(int argc, const char** argv)
{
    using namespace std;
    using namespace varoom;
    using namespace varoom::seq;

    typedef fastq_pair<fastq_reader,fastq_reader> fastq_pair_reader;
    typedef fastq_pair_reader::item_type read_pair;

    const size_t K = 25;
    const size_t N = 1000000;

    unordered_map<kmer,size_t> X;
    auto t0 = std::chrono::high_resolution_clock::now();

    vector<kmer> lhsFwd;
    vector<kmer> lhsRev;
    vector<kmer> rhsFwd;
    vector<kmer> rhsRev;

    input_file_holder_ptr inp1 = files::in(argv[1]);
    input_file_holder_ptr inp2 = files::in(argv[2]);
    fastq_reader r1(**inp1);
    fastq_reader r2(**inp2);
    for (fastq_pair_reader r1r2(r1, r2); r1r2.more(); ++r1r2)
    {
        kmers::make(std::get<1>((*r1r2).first), K, lhsFwd, lhsRev);
        kmers::make(std::get<1>((*r1r2).second), K, rhsFwd, rhsRev);

        for (auto itr = lhsFwd.begin(); itr != lhsFwd.end(); ++itr)
        {
            X[*itr] += 1;
        }
        for (auto itr = rhsFwd.begin(); itr != rhsFwd.end(); ++itr)
        {
            X[*itr] += 1;
        }
        for (auto itr = lhsRev.begin(); itr != lhsRev.end(); ++itr)
        {
            X[*itr] += 1;
        }
        for (auto itr = rhsRev.begin(); itr != rhsRev.end(); ++itr)
        {
            X[*itr] += 1;
        }
        if (X.size() >= N)
        {
            auto t1 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> d = t1 - t0;
            t0 = t1;
            std::cerr << '\r' << (double(N) / d.count());
            X.clear();
        }
    }
    std::cerr << endl;

    return 0;
}
