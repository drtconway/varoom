#include "varoom/kmers.hpp"
#include "varoom/seq/fastq.hpp"
#include "varoom/util/files.hpp"
#include "varoom/command.hpp"

#include <unordered_map>

int main(int argc, const char** argv)
{
    using namespace std;
    using namespace varoom;
    using namespace varoom::seq;

    const size_t K = 25;

    unordered_map<kmer,size_t> X;
    vector<kmer> xs;

    input_file_holder_ptr inp = files::in(argv[1]);
    for (fastq_reader r(**inp); r.more(); ++r)
    {
        const fastq_read& rd = *r;
        kmers::make(std::get<1>(rd), K, xs);
        for (auto itr = xs.begin(); itr != xs.end(); ++itr)
        {
            X[*itr] += 1;
        }
    }
    for (auto itr = X.begin(); itr != X.end(); ++itr)
    {
        cout << kmers::render(K, itr->first) << '\t' << itr->second << endl;
    }

    return 0;
}
