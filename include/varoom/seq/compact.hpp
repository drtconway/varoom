#ifndef VAROOM_SEQ_COMPACT_HPP
#define VAROOM_SEQ_COMPACT_HPP

#include <string>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/suffix_arrays.hpp>

#include <chrono>

#ifndef VAROOM_UTIL_FILES_HPP
#include "varoom/util/files.hpp"
#endif

namespace varoom
{
    namespace seq
    {
        class compact
        {
        public:
            compact(const std::string& p_basename)
            {
                sdsl::load_from_file(m_wt, p_basename);
            }

            static void make_csa(const std::string& p_src, const std::string& p_basename)
            {
                using namespace std::chrono;
                using timer = std::chrono::high_resolution_clock;

                sdsl::csa_sada<> vv;
                auto t0 = timer::now();
                sdsl::construct(vv, p_src, 1);
                auto t1 = timer::now();
                sdsl::store_to_file(vv, p_basename);
                auto t2 = timer::now();
                std::cerr << "construction: " << duration_cast<milliseconds>(t1-t0).count() << std::endl;
                std::cerr << "storing: " << duration_cast<milliseconds>(t2-t1).count() << std::endl;
                std::cerr << "extract: " << sdsl::extract(vv, 60985588, 60985608);
            }

            static void make_wt(const std::string& p_src, const std::string& p_basename)
            {
                using namespace std::chrono;
                using timer = std::chrono::high_resolution_clock;

                sdsl::wt_huff<> vv;
                auto t0 = timer::now();
                sdsl::construct(vv, p_src, 1);
                auto t1 = timer::now();
                sdsl::store_to_file(vv, p_basename);
                auto t2 = timer::now();
                std::cerr << "construction: " << duration_cast<milliseconds>(t1-t0).count() << std::endl;
                std::cerr << "storing: " << duration_cast<milliseconds>(t2-t1).count() << std::endl;
            }

            std::string operator()(size_t p_begin, size_t p_end) const
            {
                return std::string(m_wt.begin() + p_begin, m_wt.begin() + p_end);
            }

        private:
            sdsl::wt_huff<> m_wt;
        };
    }
    // namespace seq
}
// namespace varoom

#endif // VAROOM_SEQ_COMPACT_HPP
