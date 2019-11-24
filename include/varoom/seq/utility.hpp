#ifndef VAROOM_SEQ_UTILITY_HPP
#define VAROOM_SEQ_UTILITY_HPP

namespace varoom
{
    namespace seq
    {
        class utility
        {
        public:
            static void reverse_complement(const std::string& p_seq, std::string& p_res)
            {
                p_res.clear();
                for (auto itr = p_seq.rbegin(); itr != p_seq.rend(); ++itr)
                {
                    switch (*itr)
                    {
                        case 'A':
                        case 'a':
                        {
                            p_res.push_back('T');
                            break;
                        }
                        case 'C':
                        case 'c':
                        {
                            p_res.push_back('G');
                            break;
                        }
                        case 'G':
                        case 'g':
                        {
                            p_res.push_back('C');
                            break;
                        }
                        case 'T':
                        case 't':
                        {
                            p_res.push_back('A');
                            break;
                        }
                        case 'N':
                        case 'n':
                        {
                            p_res.push_back('N');
                            break;
                        }
                        default:
                        {
                            throw std::runtime_error("unimplemented base");
                        }
                    }
                }
            }
        };
    }
    // namespace seq
}
// namespace varoom

#endif // VAROOM_SEQ_UTILITY_HPP
